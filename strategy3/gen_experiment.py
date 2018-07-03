import argparse
parser = argparse.ArgumentParser(prog='gen_experiment.py', description='''
    Generate the content needed for experiment section in config of this module
    This script will grab all files in `--dir` with BAM extention
    CAUTION: please always use absolute path in `--dir`
''')
parser.add_argument('--dir')
args = parser.parse_args()

import subprocess
import glob
import ntpath
import os

# pwd = subprocess.check_output('pwd').strip()
bams = glob.glob('{dir}/*.bam'.format(dir = args.dir))
for bam in bams:
    name = os.path.splitext(ntpath.basename(bam))[0]
    out = '''  {name}:
    bam: "{bam}"'''.format(bam = bam, name = name)
    print(out)
