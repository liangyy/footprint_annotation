import argparse
parser = argparse.ArgumentParser(prog='get_score_with_fimo.py', description='''
    Compute motif score using fimo
    in strategy2
''')
parser.add_argument('--snp', help='''
    footprint_snp
''')
parser.add_argument('--seq', help='''
    seq n TAB file
''')
parser.add_argument('--ncol_of_snp_list', type = int, help = '''
    skip ncol_of_snp_list columns to read motif information in footprint_snp
''')
parser.add_argument('--out')
parser.add_argument('--motif_folder', help='''
    motif folder
''')
args = parser.parse_args()

import sys
import os
import pandas as pd
import numpy as np
import gzip

def change_allele(seq, strand, allele_start, region_start, to_char):
    off_set =  allele_start - region_start
    seq = seq.upper()
    rdic = { 'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G', 'N' : 'N' }
    seq_new = ''
    if strand == '+':
        seq_new = seq[:off_set] + to_char + seq[off_set + 1:]
    elif strand == '-':
        off_set_from_end = len(seq) - off_set - 1
        seq_new = seq[:off_set_from_end] + rdic[to_char] + seq[off_set_from_end + 1:]
    return seq_new

def split_output(df):
    i1 = []
    i2 = []
    s1 = []
    s2 = []
    for i, r in df.iterrows():
        seqid = r['seqid']
        score = float(r['score'])
        temp = seqid.split('-')
        if temp[1] == 'seq1':
            i1.append(int(temp[0]))
            s1.append(score)
        elif temp[1] == 'seq2':
            i2.append(int(temp[0]))
            s2.append(score)
    return i1, i2, s1, s2
useful_cols = [ i for i in range(5) ] + [ j for j in range(args.ncol_of_snp_list, args.ncol_of_snp_list + 5) ]
variants = pd.read_table(args.snp, sep = '\t', compression = 'gzip', usecols = useful_cols,
    header = None, names = ['chr', 'start', 'end', 'ref', 'alt', 'chr.motif', 'start.motif', 'end.motif', 'motif', 'strand'])

if variants.shape[0] == 0:
    temp = pd.DataFrame([])
    temp.to_csv(args.out,
        sep='\t',
        compression='gzip',
        header=False, 
        index=False)
    sys.exit()

sequences = pd.read_table(args.seq, sep = '\t', header = None, names = ['info', 'seq'])
variants['seq'] = sequences['seq']
variants['seq1'] = variants.apply(lambda row: change_allele(row['seq'], row['strand'], row['start'], row['start.motif'], row['ref']), axis = 1)
variants['seq2'] = variants.apply(lambda row: change_allele(row['seq'], row['strand'], row['start'], row['start.motif'], row['alt']), axis = 1)
variants['id'] = range(variants.shape[0])
variants['score1'] = np.nan
variants['score2'] = np.nan
for motif in variants.motif.unique():
    print(motif)
    print(variants)
    motif_file = '{motif_folder}/{motif_name}.meme'.format(motif_folder = args.motif_folder, motif_name = motif)
    subset = variants.loc[variants.motif == motif]
    temp_in_name = '{out}.fimo_temp.in'.format(out = args.out)
    temp_out_name = '{out}.fimo_temp.out'.format(out = args.out)
    temp_o = open(temp_in_name, 'w')
    for i, r in subset.iterrows():
        temp_o.write('> {id}-seq1\n{seq}\n'.format(id = r['id'], seq = r['seq1']))
    for i, r in subset.iterrows():
        temp_o.write('> {id}-seq2\n{seq}\n'.format(id = r['id'], seq = r['seq2']))
    temp_o.close()
    command = 'fimo --skip-matched-sequence --norc --text --thresh 1 {motif_file} {temp_in} > {temp_out}'.format(motif_file = motif_file, temp_in = temp_in_name, temp_out = temp_out_name)
    print(command)
    os.system(command)
    temp_out_df = pd.read_table(temp_out_name, sep = '\t', header = None, skiprows = 1, usecols = [0, 2, 5, 6], names = ['motif', 'seqid', 'strand', 'score'])
    id1, id2, s1, s2 = split_output(temp_out_df)
    variants.ix[id1, 'score1'] = s1
    variants.ix[id2, 'score2'] = s2
    command = 'rm {temp_in}'.format(temp_in = temp_in_name)
    os.system(command)
    command = 'rm {temp_out}'.format(temp_out = temp_out_name)
    os.system(command)
print(variants)
variants.iloc[:, [0, 1, 2, 3, 4, 9, 8, 14, 15]].to_csv(args.out, sep = '\t', header = False, index = False, compression = 'gzip')
