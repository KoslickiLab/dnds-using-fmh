import random
from more_itertools import sliced
from Bio.Seq import Seq

# divide a nucleotide sequence into its codons creating a list of codons
def get_coding_sequence_from_nucleotide_sequence(nt_sequence):
    """This functioin returns the coding sequence as a list when given a nucleotide sequence.
    nt_sequence: nucleotide sequence, preferably a protein coding gene sequence"""
    coding_sequence = list(sliced(nt_sequence,3))     #coding sequence
    codons_to_remove = ['TGA', 'TAA', 'TAG']     #list of stop codons
    filtered_codons = [codon for codon in codons if codon not in codons_to_remove]     # Remove stop codons
    return(coding_sequence)

# decision to mutate a position by desired p rate
def mutate_position_based_on_mutation_rate_p(p_mutation_rate):
    """This function returns a boolean probability of a mutation based off a given mutation rate p
    p_mutation_rate = 1 - Cfrac(a,B)**(1/k)
    Cfrac is the containment index between two sequences"""
    p = random.random()
    if p <= p_mutation_rate:
        return(True)
    else:
        return(False)

# mutate a position by a new nucleotide
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

"""
def mutated_sequence_based_on_mutation_rate_p(sequence,p_mutation_rate):
    #This function taks in a sequence to be mutated based on the mutation rate p, which is another argument of the function. 
    #The function loops through each position of the sequence,
    #decides whether the nucleotide at that position is mutated, 
    #and adds on to the newly mutated sequence, which will then be return in the end.
    mutated_sequence = ''
    for nt_position in sequence.upper():
        if mutate_position_based_on_mutation_rate_p(p_mutation_rate):
            mutate_with = random.choice(mutate_with_nucleotides(nt_position.upper()))
            mutated_sequence+=mutate_with
        else:
            mutated_sequence+=nt_position
    #print('muated seq',mutated_sequence)
    return(mutated_sequence.upper())
"""

def mutate_position(nt_position, p_mutation_rate):
    if mutate_position_based_on_mutation_rate_p(p_mutation_rate):
        mutate_with = random.choice(mutate_with_nucleotides(nt_position.upper()))
        return(mutate_with)
    else:
        return(nt_position)

def positive_selection_outcome(codon,p_mutation_rate):
    mutated_codon = ''
    if random.choice([0, 1, 2])==random.choice([0, 1, 2]): 
        for position in len(codon):
            nt_position = codon[position]
            if position == 0 or position == 1:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
            elif position == 2:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
    else:
        for position in len(codon):
            nt_position = codon[position]
            if position == 0 or position == 1:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
    return(mutated_codon)

def negative_selection_outcome(codon,p_mutation_rate):
    mutated_codon = ''
    if random.choice([0, 1, 2])==random.choice([0, 1, 2]): 
        for position in len(codon):
            nt_position = codon[position]
            if position == 0 or position == 1:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
            elif position == 2:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
    else:
        for position in len(codon):
            nt_position = codon[position]
            if position == 2:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
    return(mutated_codon)

def positive_selection_based_on_mutation_rate_p(sequence,p_mutation_rate):
    """This function taks in a sequence to be mutated based on the mutation rate p, which is another argument of the function. 
    The function loops through each position of the sequence,
    decides whether the nucleotide at that position is mutated, 
    and adds on to the newly mutated sequence, which will then be return in the end."""
    mutated_sequence = ''
    coding_sequence = get_coding_sequence_from_nucleotide_sequence(sequence)
    for codon in coding_sequence.upper():
        mutated_sequence+=positive_selection_outcome(codon, p_mutation_rate)
    #print('muated seq',mutated_sequence)
    return(mutated_sequence.upper())

def negative_selection_based_on_mutation_rate_p(sequence,p_mutation_rate):
    """This function taks in a sequence to be mutated based on the mutation rate p, which is another argument of the function. 
    The function loops through each position of the sequence,
    decides whether the nucleotide at that position is mutated, 
    and adds on to the newly mutated sequence, which will then be return in the end."""
    mutated_sequence = ''
    coding_sequence = get_coding_sequence_from_nucleotide_sequence(sequence)
    for codon in coding_sequence.upper():
        mutated_sequence+=negative_selection_outcome(codon, p_mutation_rate)
    #print('muated seq',mutated_sequence)
    return(mutated_sequence.upper())

