import argparse
parser = argparse.ArgumentParser(prog='get_motif_scores.py', description='''
    For each snp listed in the input bed file, extract the corresponding genome
    sequence with reference and alternative allele. Calcalate the motif score in
    the corresponding strand direction for both reference allele and alternative
    one along with the estimated prior probability of binding by interpolation
    (get the points from motif file)
''')
parser.add_argument('--motif_folder', help='''
    Motif folder that contains the motif files which should match with the input
    footprint SNP file
''')
parser.add_argument('--footprint_snp', help='''
    Footprint SNP file obtained by intersectin SNP list with files in footprint
    region files downloaded here, http://genome.grid.wayne.edu/centisnps/bytissue/
''')
parser.add_argument('--footprint_seq', help='''
    Sequence TAB file containing sequence extracted from reference genome
''')
parser.add_argument('--ncol', type=int, help='''
    Number of columns in original input SNP list file
''')
parser.add_argument('--out')
args = parser.parse_args()

import sys
import os
if '../scripts/' not in sys.path:
    sys.path.insert(0, '../scripts/')
import footprint_lib
import gzip

f = open(args.footprint_seq,'r')
file_content = f.readlines()
f.close()
ref_seqs = []
for i in file_content:
    i = i.strip().upper()
    ref_seqs.append(i.split('\t')[-1])

o = gzip.open(args.out, 'wb')
o.write('\t'.join(['SNP.ID', 'LLR.Ref', 'LLR.Alt', 'Prior.Ref', 'Prior.Alt']) + '\n')
with gzip.open(args.footprint_snp,'rb') as f:
    counter = 0
    for line in f:
        snp, region = footprint_lib.read_bed_line(line.decode(), args.ncol)
        ref, alt = footprint_lib.get_seq(snp, region, ref_seqs[counter])
        counter += 1
        if ref is None:
            continue
        llrs = []
        priors = []
        for i in (ref, alt):
            motif = footprint_lib.get_motif(args.motif_folder, region.motif)
            llr = footprint_lib.motif_score(i, motif)
            prior = footprint_lib.bind_prior(llr, motif)
            llrs.append(llr)
            priors.append(prior)
        o.write('{idx}\t{llr_ref}\t{llr_alt}\t{prior_ref}\t{prior_alt}'.format(idx=snp.idx,
                                                                                llr_ref=llrs[0],
                                                                                llr_alt=llrs[1],
                                                                                prior_ref=priors[0],
                                                                                prior_alt=priors[1]) + '\n')
o.close()
