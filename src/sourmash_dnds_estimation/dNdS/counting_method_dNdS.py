
"""counting methods used in dN/dS estimation for codon substitution rates"""
import math

def get_coding_sequence_from_nucleotide_sequence(nt_sequence):
    """This functioin returns the coding sequence as a list when given a nucleotide sequence.
    nt_sequence: nucleotide sequence, preferably a protein coding gene sequence"""
    return(list(sliced(nt_seq[frame-1:].upper(),3)))

def total_synonymous_sites(nt_sequence):
    """Count the total number of possible syonymous mutations.
    Change nt_sequence string to a codon sequence list.
    Loop through each codon sequence and then loop through each 
    of the three nt within the codon to record syonynmous sites
    nt_sequence: a nucleotide sequence that is a string"""
    return('Total number of synonymous sites of two sequences')

def total_nonsynonymous_sites(nt_sequence):
    """Count the total number of possible nonsyonymous mutations
    Change nt_sequence string to a codon sequence list.
    Loop through each codon sequence and then loop through each 
    of the three nt within the codon to record nonsyonynmous sites
    nt_sequence: a nucleotide sequence that is a string"""
    return('Total number of nonsynonymous sites of two sequences')

def total_synonymous_diff(nt_sequence_1,nt_sequence_2):
    """Count the total number of possible syonymous differencess. 
    This is a pairwise comparison. Change both nt sequences into a codon list.
    Loop through one codon sequence and have a nested loop for the second codon sequence.
    Record syonymous differences.
    nt_sequence_1: a nucleotide sequence that is a string
    nt_sequence_2: a nucleotide sequence that is a string"""
    return('Total number of synonymous differences between two sequences')

def total_nonsynonymous_diff(nt_sequence_1,nt_sequence_2):
    """Count the total number of possible nonsyonymous differencess. 
    This is a pairwise comparison. Change both nt sequences into a codon list.
    Loop through one codon sequence and have a nested loop for the second codon sequence.
    Record nonsyonymous differences.
    nt_sequence_1: a nucleotide sequence that is a string
    nt_sequence_2: a nucleotide sequence that is a string"""
    return('Total number of nonsynonymous differences between two sequences')

def synonymous_diff_proportions(nt_sequence_1,nt_sequence_2):
    """Obtains synonymous difference proportions.
    nt_sequence_1: a nucleotide sequence that is a string
    nt_sequence_2: a nucleotide sequence that is a string"""
    pS = total_synonymous_diff(nt_sequence_1,nt_sequence_2)/total_synonymous_sites(nt_sequence)
    return(pS)

def nonsynonymous_diff_proportions(nt_sequence_1,nt_sequence_2):
    """Obtains nonsynonymous difference proportions.
    nt_sequence_1: a nucleotide sequence that is a string
    nt_sequence_2: a nucleotide sequence that is a string"""
    pN = total_nonsynonymous_diff(nt_sequence_1,nt_sequence_2)/total_nonsynonymous_sites(nt_sequence)
    return(pN)

#NG86 is one of the counting methods
def JC69_correction(diiference_proportion_value):
    """Reports the JC69 correction and produces either a corrected dS or dN"""
    JC69=(-3/4)*(math.log(1-(4/3)*(diiference_proportion_value)))
    return(JC69)

def NG86(nt_sequence_1,nt_sequence_2):
    dN=JC69_correction(nonsynonymous_diff_proportions(nt_sequence_1,nt_sequence_2))
    dS=JC69_correction(synonymous_diff_proportions(nt_sequence_1,nt_sequence_2))
    return(dN/dS)

#def K80_correction():

#def LWL85():

#def LPB93():

