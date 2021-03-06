import argparse
parser = argparse.ArgumentParser(prog='scan_motif_in_region.py', description='''
    This script takes the regions of interest and motif list. It computes the log
    likelihood ratio under certain background base distribution.
''')
parser.add_argument('--fasta', help='''
    Sequence file with name telling the genomic location of the region
''')
parser.add_argument('--motif_dir', help='''
    Motif directory
''')
parser.add_argument('--motif_name', help='''
    Motif name in motif directory (naming convention: [name].pwm)
''')
parser.add_argument('--output', help='''
    Name of output file (it is in bed format with the last two columns telling
    the motif name and the llr score)
''')
# parser.add_argument('--method', help='''
#     The method to call 'active region' namely motif binding site.
# ''')
parser.add_argument('--threshold', type=float, help='''
    The threshold used to call 'active region'. It should be used along with
    --method
''')

args = parser.parse_args()

import sys, re, gzip
if '../scripts/' not in sys.path:
    sys.path.append('../scripts/')
import footprint_lib
import numpy as np


motif = footprint_lib.get_motif(args.motif_dir, args.motif_name)
motif_length = motif.pwm.shape[0]
out = gzip.open(args.output, 'w')
with gzip.open(args.fasta, 'r') as f:
    for i in f:
        i = i.decode()
        i = i.strip()
        if i[0] == '>':
            print(i)
            chrm = re.search('>(chr[0-9XY]+):', i).group(1)
            start = int(re.search('>chr[0-9XY]+:([0-9]+)-', i).group(1))
            continue
        else:
            seq = i
            passed_regions = footprint_lib.scan_region(seq, motif, np.log2(args.threshold))
        for r in passed_regions:
            pos, score, strand = r
            rstart = str(start + pos)
            rend = str(start + pos + motif_length)
            new_line = '\t'.join([chrm, rstart, rend, strand, motif.name, str(np.round(score, decimals = 4))]) + '\n'
            out.write(new_line.encode('ascii')) 
out.close()
