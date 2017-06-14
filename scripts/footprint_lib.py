import os, sys
import numpy as np

class SNP:
    def __init__(self, chrm, start, end, ref, alt, idx):
        self.chr = chrm
        self.start = start
        self.end = end
        self.ref = ref
        self.alt = alt
        self.idx = idx

class Footprint:
    def __init__(self, chrm, start, end, motif, strand, bind_prior):
        self.chr = chrm
        self.start = start
        self.end = end
        self.motif = motif
        self.strand =  strand
        self.bind_prior = bind_prior

class Motif:
    def __init__(self, pwm, priors, llrs):
        self.pwm = pwm
        self.priors = priors
        self.llrs = llrs

def read_bed_line(line, ncol):
    line = line.split('\t')
    snp = SNP(line[0], int(line[1]), int(line[2]), line[3], line[4], line[5])
    line = line[ncol : ]
    footprint = Footprint(line[0], int(line[1]), int(line[2]), line[3], line[4], float(line[5]))
    return snp, footprint

def get_motif(dirname, motif_name):
    f = open(os.sep.join([dirname, '{motif_name}.pwm'.format(motif_name=motif_name)]), 'r')
    lines = f.readlines()
    priors = lines[1].strip().split('\t')[-9:]
    priors = [ float(i) for i in priors ]
    llrs = lines[3].strip().split('\t')[-9:]
    llrs = [ float(i) for i in llrs ]
    lines = lines[5:]
    pwm = []
    for line in lines: # A, C, G, T in order
        pwm.append([ float(i) for i in line.split('\t') ])
    pwm = np.array(pwm)
    motif = Motif(pwm, priors, llrs)
    return motif

def motif_score(seq, motif): # seq is in character
    # digit_seq = _to_digit(seq)
    llr_mat = np.log2(motif.pwm + 1e-6) - np.log2(np.ones(motif.pwm.shape) / 4)
    llr = (seq * llr_mat).sum()
    return llr

def bind_prior(llr, motif):
    xp = motif.llrs
    fp = motif.priors
    prior = np.interp(llr, xp, fp)
    return prior

def _to_digit(seq):
    digit = np.zeros((len(seq), 4))
    dic = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }
    for i in range(len(seq)):
        digit[i, dic[seq[i]]] = 1
    return digit

def get_seq(snp, region, seq):
    # seq = _get_from_fasta(region.chr, region.start, snp.end, genome)
    ref = seq
    pos = snp.start - region.start
    # print(pos, snp.start, region.start)
    if seq[pos] != snp.ref:
         print('Ref in SNP is {ref_snp} does not match Ref in fasta {ref_fa}. Skip!'.format(ref_snp=snp.ref, ref_fa=seq[pos]), file=sys.stderr)
         return None, None
    alt = seq[:pos] + snp.alt + seq[pos+1:]
    ref = _to_digit(ref)
    alt = _to_digit(alt)
    if region.strand == '-':
        ref = ref[::-1,::-1]
        alt = alt[::-1,::-1]
    return ref, alt
