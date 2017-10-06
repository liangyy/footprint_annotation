import argparse
parser = argparse.ArgumentParser(prog='pwm2meme.py', description='''
    Convert pwm motif into meme format
''')
parser.add_argument('--motif_dir_in_pwm', help='''
    Motif directory of input motif
''')
parser.add_argument('--motif_name_in_pwm', help='''
    Motif name of input motif
''')
parser.add_argument('--output_in_meme')
args = parser.parse_args()

import sys, re, gzip
if '../scripts/' not in sys.path:
    sys.path.append('../scripts/')
import footprint_lib
import numpy as np

motif = footprint_lib.get_motif(args.motif_dir_in_pwm, args.motif_name_in_pwm)
out = open(args.output_in_meme, 'w')
width = motif.pwm.shape[0]
header = '''MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies
A 0.25 C 0.25 G 0.25 T 0.25

MOTIF {motif_name}
letter-probability matrix: alength= 4 w= {length}'''.format(
                                                            length = width,
                                                            motif_name = args.motif_name_in_pwm)
out.write(header)
for i in range(width):
    for j in motif.pwm[i]:
        out.write(str(j))
        out.write(' ')
    out.write('\n')
