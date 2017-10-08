import argparse
parser = argparse.ArgumentParser(prog='thresholding_by_motif_prior.py', description='''
    Filter out score that are too small by motif prior prob
''')
parser.add_argument('--input')
parser.add_argument('--motif_dir')
parser.add_argument('--motif_name')
parser.add_argument('--output')
args = parser.parse_args()

import sys, os
if '../scripts/' not in sys.path:
    sys.path.append('../scripts/')

import footprint_lib

motif = footprint_lib.get_motif(args.motif_dir, args.motif_name)
for i in range(len(motif.priors)):
    if motif.priors[i] == 0.1:
        llr = motif.llrs[i]
cmd = '''
    zcat {input} | awk -F"\\t" '$5 > {llr}' | gzip > {output}
'''.format(input = args.input, llr = llr, output = args.output)
print(cmd)
os.system(cmd)
