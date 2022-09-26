from helperfuncs import *

def determine_approx_dnds_no_fmh(dna_seq1, dna_seq2, k):
    nt_containment = sequences_to_containment_perfect(dna_seq1, dna_seq2, 3*k)
    T = translate_dna_to_aa
    aa_containemnt = sequences_to_containment_perfect( T(dna_seq1), T(dna_seq2), k )
    p_nt = containment_to_mut_rate(nt_containment, 3*k)
    p_aa = containment_to_mut_rate(aa_containemnt, k)

    return p_aa/(1.0 - p_aa - (1 - p_nt)**3)


# given two related dna sequences, determine dn/ds
def determine_approx_dnds_using_fmh(dna_seq1, dna_seq2, k, scale_factor):
    aa_k = k
    nt_k = 3*k
    nt_containment = sequences_to_containment_using_fmh(dna_seq1, dna_seq2, nt_k, scale_factor)
    aa_containemnt = sequences_to_containment_using_fmh( translate_dna_to_aa(dna_seq1), translate_dna_to_aa(dna_seq2), aa_k, scale_factor)
    #print(nt_containment, aa_containemnt)
    p_nt = containment_to_mut_rate(nt_containment, nt_k)
    p_aa = containment_to_mut_rate(aa_containemnt, aa_k)

    try:
        dnds = p_aa / ( 1.0 - p_aa - (1.0 - p_nt)**3 )
    except:
        if p_nt == p_aa:
            if p_nt == 0.0:
                print('No mutations! Identical.')
                return -1
            else:
                dnds = float('Infinity')
        else:
            print('Some unknown error occurred!')
            return -2
    return dnds

def determine_correct_dnds(dna_seq1, dna_seq2):
    num_total_mutations = 0
    for i in range(len(dna_seq1)):
        if dna_seq1[i] != dna_seq2[i]:
            num_total_mutations += 1
    aa_seq1 = translate_dna_to_aa(dna_seq1)
    aa_seq2 = translate_dna_to_aa(dna_seq2)
    num_non_syn_mutations = 0
    for i in range(len(aa_seq1)):
        if aa_seq1[i] != aa_seq2[i]:
            num_non_syn_mutations += 1
    #print(num_total_mutations, num_non_syn_mutations)
    if num_total_mutations == num_non_syn_mutations:
        if num_total_mutations == 0:
            print('No mutations! Identical.')
            return -1
        else:
            # all mutations are nonsynonymous
            return float('Infinity')

    return 1.0 * num_non_syn_mutations / (num_total_mutations-num_non_syn_mutations)
