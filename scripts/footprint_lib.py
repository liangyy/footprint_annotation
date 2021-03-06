import os, sys
import numpy as np
from scipy import stats

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
    def __init__(self, pwm, priors, llrs, motif_name):
        self.pwm = pwm
        self.priors = priors
        self.llrs = llrs
        self.name = motif_name

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
    llrs = lines[2].strip().split('\t')[-9:]
    llrs = [ float(i) for i in llrs ]
    lines = lines[5:]
    pwm = []
    for line in lines: # A, C, G, T in order
        pwm.append([ float(i) for i in line.split('\t') ])
    pwm = np.array(pwm)
    motif = Motif(pwm, np.array(priors), np.array(llrs), motif_name)
    return motif

def motif_score(seq, motif): # seq is in character
    # digit_seq = _to_digit(seq)
    llr_mat = np.log2(motif.pwm + 1e-6) - np.log2(np.ones(motif.pwm.shape) / 4)
    llr = (seq * llr_mat).sum()
    return llr

def bind_prior(llr, motif):
    y = np.log(motif.priors) - np.log(1 - motif.priors)
    x = motif.llrs
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    r_square = r_value ** 2
    if r_square < 0.99:
        print('R^2 = {r} is less than 0.99. Skip this motif {name}'.format(r=r_square, name=motif.name), file=sys.stderr)
        return None
    log_ratio_prior = llr * slope + intercept
    return log_ratio_prior

def _to_digit(seq):
    seq = seq.upper()
    digit = np.zeros((len(seq), 4))
    dic = { 'A': 0, 'C': 1, 'G': 2, 'T': 3 }
    for i in range(len(seq)):
        if seq[i] in dic.keys():
            digit[i, dic[seq[i]]] = 1
        else:
            digit[i, :] = 0.25
    return digit

def get_seq(snp, region, seq, check_ref):
    # seq = _get_from_fasta(region.chr, region.start, snp.end, genome)
    ref = seq
    pos = snp.start - region.start
    # print(pos, snp.start, region.start)
    if seq[pos] != snp.ref and check_ref == '1':
         print('Ref in SNP is {ref_snp} does not match Ref in fasta {ref_fa}. Skip!'.format(ref_snp=snp.ref, ref_fa=seq[pos]), file=sys.stderr)
         return None, None, None
    ref = seq[:pos] + snp.ref + seq[pos+1:]
    alt = seq[:pos] + snp.alt + seq[pos+1:]
    ref = _to_digit(ref)
    alt = _to_digit(alt)
    if region.strand == '-':
        ref = ref[::-1,::-1]
        alt = alt[::-1,::-1]
        pos = region.end - region.start - pos
    else:
        pos = pos + 1
    return ref, alt, pos

def scan_region(seq, motif, threshold):
    seq_digit = _to_digit(seq)
    seq_length = len(seq)
    motif_length = motif.pwm.shape[0]
    out = []
    for i in range(seq_length - motif_length + 1):
        subseq_digit = seq_digit[i : i + motif_length, :]
        subseqrev_digit = subseq_digit[::-1, ::-1]
        forward_score = motif_score(subseq_digit, motif)
        reverse_score = motif_score(subseqrev_digit, motif)
        if forward_score > threshold:
            out.append([i, forward_score, '+'])
        if reverse_score > threshold:
            out.append([i, reverse_score, '-'])
    return out
