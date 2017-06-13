import argparse
parser = argparse.ArgumentParser(prog='get_motif_scores.py', description='''
    For each snp listed in the input bed file, extract the corresponding genome
    sequence with reference and alternative allele. Calcalate the motif score in
    the corresponding strand direction for both reference allele and alternative
    one along with the estimated prior probability of binding by interpolation
    (get the points from motif file)
    ''')
parser.add_argument('--genome', help='''
    Genome fasta file that should extract sequence using coordinates
''')
parser.add_argument('--motif_folder', help='''
    Motif folder that contains the motif files which should match with the input
    footprint SNP file
''')
parser.add_argument('--footprint_snp', help='''
    Footprint SNP file obtained by intersectin SNP list with files in footprint
    region files downloaded here, http://genome.grid.wayne.edu/centisnps/bytissue/
''')
parser.add_argument('--ncol', type=int, help='''
    Number of columns in original input SNP file
''')
args = parser.parse_args()

import sys
import os
if '../scipts/' not in sys.path:
    sys.path.insert(0, 'scripts/')
import footprint_lib

import gzip

with gzip.open(args.footprint_snp,'r') as f:
    for line in f:
        snp, region = footprint_lib.read_bed_line(line, args.ncol)
        ref, alt = footprint_lib.get_seq(snp, region, args.genome)
        if ref is None:
            continue
        llrs = []
        priors = []
        for i in (ref, alt):
            motif = footprint_lib.get_motif(region.motif, args.motif_folder)
            llr = footprint_lib.motif_score(i, motif)
            prior = footprint_lib.bind_prior(llr, motif)
            llrs.append(llr)
            priors.append(prior)
        print('{idx}\t{score_ref}\t{score_alt}\t{prior_ref}\t{prior_alt}'.format(idx=snp.idx,
                                                                                llr_ref=llrs[0],
                                                                                llr_alt=llrs[1],
                                                                                prior_ref=priors[0],
                                                                                prior_alt=priors[1]))
