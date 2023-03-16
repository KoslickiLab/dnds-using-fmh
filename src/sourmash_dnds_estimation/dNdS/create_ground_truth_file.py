import random
from more_itertools import sliced
from Bio.Seq import Seq

def codon_table():
    """When invoked, function returns the codon table.
    A codon table is a dictionary where the key is the codon and the values are the amino acid"""
    codontab = {
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'AGC': 'S', 'AGT': 'S',    # Serine
        'TTC': 'F', 'TTT': 'F',                                                    # Phenilalanine
        'TTA': 'L', 'TTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',    # Leucine
        'TAC': 'Y', 'TAT': 'Y',                                                    # Tirosine
        'TAA': '*', 'TAG': '*', 'TGA': '*',                                        # Stop
        'TGC': 'C', 'TGT': 'C',                                                    # Cisteine
        'TGG': 'W',                                                                # Tryptophane
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',                            # Proline
        'CAC': 'H', 'CAT': 'H',                                                    # Histidine
        'CAA': 'Q', 'CAG': 'Q',                                                    # Glutamine
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',                            # Arginine
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I',                                        # Isoleucine
        'ATG': 'M',                                                                # Methionine
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',                            # Threonine
        'AAC': 'N', 'AAT': 'N',                                                    # Asparagine
        'AAA': 'K', 'AAG': 'K',                                                    # Lysine
        'AGA': 'R', 'AGG': 'R',                                                    # Arginine
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',                            # Valine
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',                            # Alanine
        'GAC': 'D', 'GAT': 'D',                                                    # Aspartic Acid
        'GAA': 'E', 'GAG': 'E',                                                    # Glutamic Acid
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'                             # Glycine
    }
    return(codontab)

"""
def codon_substitution():
    codon_sub = {
        'TCA': {'synonymous':1,'nonsynonymous':2},    # Serine
        'TCC': {'synonymous':1,'nonsynonymous':2}, 
        'TCG': {'synonymous':1,'nonsynonymous':2}, 
        'TCT': {'synonymous':1,'nonsynonymous':2}, 
        'AGC': {}, 
        'AGT': {},
        'TTC': {'synonymous':,'nonsynonymous':},    # Phenilalanine 
        'TTT': {'synonymous':,'nonsynonymous':},
        'TTA': {'synonymous':,'nonsynonymous':},                                   # Leucine
        'TTG': {'synonymous':,'nonsynonymous':}, 
        'CTA': {'synonymous':,'nonsynonymous':}, 
        'CTC': {'synonymous':,'nonsynonymous':}, 
        'CTG': {'synonymous':,'nonsynonymous':}, 
        'CTT': {'synonymous':,'nonsynonymous':},    
        'TAC': {'synonymous':,'nonsynonymous':}, # Tirosine
        'TAT': {'synonymous':,'nonsynonymous':},                                                    
        'TGC': {'synonymous':,'nonsynonymous':}, # Cisteine
        'TGT': {'synonymous':,'nonsynonymous':},                                                    
        'TGG': {'synonymous':,'nonsynonymous':},                                                                # Tryptophane
        'CCA': {'synonymous':,'nonsynonymous':}, # Proline
        'CCC': {'synonymous':,'nonsynonymous':}, 
        'CCG': {'synonymous':,'nonsynonymous':}, 
        'CCT': {'synonymous':,'nonsynonymous':},                            
        'CAC': {'synonymous':,'nonsynonymous':}, # Histidine
        'CAT': {'synonymous':,'nonsynonymous':},                                                    
        'CAA': {'synonymous':,'nonsynonymous':}, # Glutamine
        'CAG': {'synonymous':,'nonsynonymous':},                                                   
        'CGA': {'synonymous':,'nonsynonymous':}, # Arginine
        'CGC': {'synonymous':,'nonsynonymous':}, 
        'CGG': {'synonymous':,'nonsynonymous':}, 
        'CGT': {'synonymous':,'nonsynonymous':},                            
        'ATA': {'synonymous':,'nonsynonymous':}, # Isoleucine
        'ATC': {'synonymous':,'nonsynonymous':}, 
        'ATT': {'synonymous':,'nonsynonymous':},                                       
        'ATG': {'synonymous':,'nonsynonymous':},                                                                # Methionine
        'ACA': {'synonymous':,'nonsynonymous':}, # Threonine
        'ACC': {'synonymous':,'nonsynonymous':}, 
        'ACG': {'synonymous':,'nonsynonymous':}, 
        'ACT': {'synonymous':,'nonsynonymous':},                            
        'AAC': {'synonymous':,'nonsynonymous':}, # Asparagine
        'AAT': {'synonymous':,'nonsynonymous':},                                                    
        'AAA': {'synonymous':,'nonsynonymous':}, # Lysine
        'AAG': {'synonymous':,'nonsynonymous':},                                                    
        'AGA': {'synonymous':,'nonsynonymous':}, # Arginine
        'AGG': {'synonymous':,'nonsynonymous':},                                                    
        'GTA': {'synonymous':,'nonsynonymous':}, # Valine
        'GTC': {'synonymous':,'nonsynonymous':}, 
        'GTG': {'synonymous':,'nonsynonymous':}, 
        'GTT': {'synonymous':,'nonsynonymous':},                            
        'GCA': {'synonymous':,'nonsynonymous':}, # Alanine
        'GCC': {'synonymous':,'nonsynonymous':}, 
        'GCG': {'synonymous':,'nonsynonymous':}, 
        'GCT': {'synonymous':,'nonsynonymous':},                            
        'GAC': {'synonymous':,'nonsynonymous':}, # Aspartic Acid
        'GAT': {'synonymous':,'nonsynonymous':},                                                    
        'GAA': {'synonymous':,'nonsynonymous':}, # Glutamic Acid
        'GAG': {'synonymous':,'nonsynonymous':},                                                    
        'GGA': {'synonymous':,'nonsynonymous':}, # Glycine
        'GGC': {'synonymous':,'nonsynonymous':}, 
        'GGG': {'synonymous':,'nonsynonymous':}, 
        'GGT': {'synonymous':,'nonsynonymous':}                             
    }
    return(codon_sub)
"""

def mutate_position_based_on_mutation_rate_p(p_mutation_rate):
    """This function returns a boolean probability of a mutation based off a given mutation rate p
    p_mutation_rate = 1 - Cfrac(a,B)**(1/k)
    Cfrac is the containment index between two sequences"""
    p = random.random()
    if p <= p_mutation_rate:
        return(True)
    elif p >= 1-p_mutation_rate:
        return(False)

def mutate_with_nucleotides(mutate_nucleotide):
    """If the probability is that a mutation occurs, 
    then this function makes sure that the random choice of a nucloetide mutation is not the same nucleotide of the position being changed.
    This function removes the nucleotide that will be changed from the nucleotide mutation choices
    and returns a list of three nucleotides that does not include the nucleotide to be mutated.
    The function takes in a nucleotide. 
    nucleotide: the nucleotide to be mutated"""
    nucleotides = ['A','G','C','T']
    nucleotides.remove(mutate_nucleotide)
    return(nucleotides)

def mutated_sequence_based_on_mutation_rate_p(sequence,p_mutation_rate):
    """This function taks in a sequence to be mutated based on the mutation rate p, which is another argument of the function. 
    The function loops through each position of the sequence,
    decides whether the nucleotide at that position is mutated, 
    and adds on to the newly mutated sequence, which will then be return in the end."""
    mutated_sequence = ''
    for nt_position in sequence.upper():
        if mutate_position_based_on_mutation_rate_p(p_mutation_rate):
            mutate_with = random.choice(mutate_with_nucleotides(nt_position.upper()))
            mutated_sequence+=mutate_with
        else:
            mutated_sequence+=nt_position
    #print('muated seq',mutated_sequence)
    return(mutated_sequence.upper())

def get_coding_sequence_from_nucleotide_sequence(nt_sequence):
    """This functioin returns the coding sequence as a list when given a nucleotide sequence.
    nt_sequence: nucleotide sequence, preferably a protein coding gene sequence"""
    return(list(sliced(nt_sequence,3)))

def translate_coding_sequence(cds_seq):
    """Function translates open reading frames using the coding sequence.
    cds_seq: a list of codons for a given sequemce"""
    codontab=codon_table()
    translation=''
    for codon in cds_seq:
        if codon in codontab:
            translation+=codontab[codon]
    return(translation)

def total_nucleotide_mutations(nt_sequence_1,nt_sequence_2):
    """Returns total number of mutations *both synonymous and nonsynonymous mutations
    total_muts records the total differences between two sequences.
    We are assuming that they should be the same length.
    nt_sequence_1: a nucleotide sequence that is a string
    nt_sequence_2: a nucleotide sequence that is a string"""
    total_nt_muts_count=0
    for i in range(len(nt_sequence_1)):
        seq_1_position_temp = nt_sequence_1[i]
        seq_2_position_temp = nt_sequence_2[i]
        if seq_1_position_temp != seq_2_position_temp:
            total_nt_muts_count+=1
    return(total_nt_muts_count)

def total_aa_differences(nt_sequence_1,nt_sequence_2): 
    """Returns total amino acid changes or nonsynonymous differences between two protein sequences
    We are assuming that they should be the same length.
    nt_sequence_1: a nucleotide sequence that is a string
    nt_sequence_2: a nucleotide sequence that is a string"""
    total_aa_muts_count=0
    #cds_1 = get_coding_sequence_from_nucleotide_sequence(nt_sequence_1)
    #cds_2 = get_coding_sequence_from_nucleotide_sequence(nt_sequence_2)
    #aa_seq_1 = translate_coding_sequence(cds_1)
    #aa_seq_2 = translate_coding_sequence(cds_2)
    aa_seq_1 = Seq(nt_sequence_1).translate()
    aa_seq_2 = Seq(nt_sequence_2).translate()
    for i in range(len(aa_seq_1)):
        seq_1_position_temp = aa_seq_1[i]
        seq_2_position_temp = aa_seq_2[i]
        if seq_1_position_temp != seq_2_position_temp:
                total_aa_muts_count+=1
    return(total_aa_muts_count)

def total_synonymous_mutations(total_nt_mutations,total_nonsyn_mutations):
    """Return total synonymous mutations
    total_nt_mutations: total nucleotide mutations between two sequences, obtain using total_nucleotide_mutations(nt_sequence_1,nt_sequence_2)
    total_nonsyn_mutations: total nucleotide mutations between two sequences, obtain using total_aa_differences(nt_sequence_1,nt_sequence_2)
    """
    total_syn_mutations = total_nt_mutations - total_nonsyn_mutations
    return(total_syn_mutations)

def koslicki_dnds(total_nonsyn_mutations,total_syn_mutations):
    """Return dN/dS estimation where dN is the 
    total number of nonsynonymous mutations divided by protein sequence length
    and dS is the total number of synonymous mutations divided by the protein sequence length,
    then calculates the ratio of dN/dS
    total_nonsyn_mutations: total nucleotide mutations between two sequences, obtain using total_aa_differences(nt_sequence_1,nt_sequence_2)
    total_syn_mutations: total nucleotide mutations between two sequences, obtain using total_synonymous_mutations(total_nt_mutations,total_nonsyn_mutations)
    protein_length: lenght of the protein sequence (assume both sequences are the same length)
    """
    dN = float(total_nonsyn_mutations)
    dS = float(total_syn_mutations)
    if dS != float(0):
        dnds = dN/dS
    else:
        dnds=None
    return(dnds)
