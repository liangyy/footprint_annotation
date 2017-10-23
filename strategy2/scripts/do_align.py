import argparse
parser = argparse.ArgumentParser(prog='do_align.py', description='''
    Merge count by region
''')
parser.add_argument('--cmd')
parser.add_argument('--reference')
parser.add_argument('--threads')
parser.add_argument('--o')
parser.add_argument('--i1')
parser.add_argument('--i2')
args = parser.parse_args()

import os, re

cmd = re.sub('\\\\', '', args.cmd)
cmd = cmd.format(ref = args.reference, threads = args.threads, sam = args.o, fastq1 = args.i1, fastq2 = args.i2)
os.system(cmd)
