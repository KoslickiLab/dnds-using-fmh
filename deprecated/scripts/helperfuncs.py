from scipy.stats import norm as scipy_norm
from Bio.Seq import Seq
import mmh3
from fracminhash import FracMinHash
import random

try:
    from mpmath import mp as mpmath,mpf
    mpmath.dps = 50
except ModuleNotFoundError:
    mpf = lambda v:float(v)

# z_alpha
def z_alpha(a):
    return scipy_norm.ppf(1.0 - a/2.0)

# mutation rate to q = 1 - (1-p)^k
def mu_rate_to_q(k, p):
    p = mpf(p)
    q = 1-(1-p)**k
    return float(q)

# expectation of Nmut
def exp_n_mutated(L, k, p):
    q = mu_rate_to_q(k, p)
    return L*q

# variance of nMut
def var_n_mutated(L, k, p):
	if (p == 0): return 0.0
	p = float(p)
	q = mu_rate_to_q(k, p)
	varN = L*(1-q)*(q*(2*k+(2/p)-1)-2*k) \
	     + k*(k-1)*(1-q)**2 \
	     + (2*(1-q)/(p**2))*((1+(k-1)*(1-q))*p-q)
	assert (varN>=0.0)
	return float(varN)

# given two strings, code them to mutated bit string
def string_to_mutated(string_1, string_2):
    mutated = []
    for i in range(len(string_1)):
        if i >= len(string_2):
            mutated.append(1)
            continue
        if string_1[i] == string_2[i]:
            mutated.append(0)
        else:
            mutated.append(1)
    return mutated

# given mutated bit string, determine containment
def mutated_bits_to_containment(mutated, k):
    if k >= len(mutated):
        if sum(mutated) == 0:
            return 1.0
        else:
            return 0.0
    num_kmers = len(mutated) - k + 1
    num_mutated = 0
    for i in range(len(mutated) - k + 1):
        if sum(mutated[i : i+k]) > 0:
            num_mutated += 1
    return 1.0 - num_mutated/num_kmers

# given dna seq, translate to AA sequence
def translate_dna_to_aa(genome_sequence):
    dna_seq = Seq(genome_sequence)
    aa_seq = dna_seq.translate()
    return str(aa_seq)

# given containment_index, return mutation rate
def containment_to_mut_rate(containment, k):
    return 1.0 - containment ** (1.0/k)

def get_hash_from_kmer(kmer, seed=0):
	hash_value = mmh3.hash64(str(kmer), seed=seed)[0]
	if hash_value < 0:
		hash_value += 2**64
	return hash_value

def sequence_to_fmh_sketch(seq, k, scale_factor):
    fmh = FracMinHash(scale_factor, 2**64)
    for i in range(len(seq) - k + 1):
        hash_value = get_hash_from_kmer( seq[i: i+k] )
        fmh.add_value(hash_value)
    return fmh

# given two sequences, determine the mutation rate
def sequences_to_containment_perfect(seq1, seq2, k):
    mutated_bits = string_to_mutated(seq1, seq2)
    containment_index = mutated_bits_to_containment(mutated_bits, k)
    #mut_rate = containment_to_mut_rate(containment_index, k)
    return containment_index

def sequences_to_mut_rate_perfect(seq1, seq2, k):
    containment = sequences_to_containment_perfect(seq1, seq2, k)
    return 1.0 - containment**(1.0/k)

# given two sequences, determine the mutation rate
def sequences_to_containment_using_fmh(seq1, seq2, k, scale_factor):
    fmh1 = sequence_to_fmh_sketch(seq1, k, scale_factor)
    fmh2 = sequence_to_fmh_sketch(seq2, k, scale_factor)
    return fmh1.get_containment(fmh2)

def sequences_to_mutrate_using_fmh(seq1, seq2, k, scale_factor):
    c = sequences_to_containment_using_fmh(seq1, seq2, k, scale_factor)
    return 1.0 - c ** (1.0 / k)

def mutate_string(in_str, mut_rate):
    out_str = list(in_str)
    other_bases = {
        'A': ['T', 'C', 'G'],
        'C': ['T', 'A', 'G'],
        'G': ['T', 'C', 'A'],
        'T': ['A', 'C', 'G']
    }
    for i in range(len(out_str)):
        if random.uniform(0,1) < mut_rate:
            out_str[i] = random.choice( other_bases[out_str[i]] )
        else:
            continue
    out_str = "".join(out_str)
    return out_str
