import random

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

def mutate_position_based_on_mutation_rate_p(p_mutation_rate):
    """This function returns a boolean probability of a mutation based off a given mutation rate p
    p_mutation_rate = 1 - Cfrac(a,B)**(1/k)
    Cfrac is the containment index between two sequences"""
    p = random.random()
    if p <= p_mutation_rate: 
        return(True)
    elif p >= 1-p_mutation_rate:
        return(False)

def mutate_with_nucleotides(nucleotide):
    """If the probability is that a mutation occurs, 
    then this function makes sure that the random choice of a nucloetide mutation is not the same nucleotide of the position being changed.
    This function removes the nucleotide that will be changed from the nucleotide mutation choices
    and returns a list of three nucleotides that does not include the nucleotide to be mutated.
    The function takes in a nucleotide. 
    nucleotide: the nucleotide to be mutated"""
    nucleotides = ['A','G','C','T']
    return(nucleotides.remove(nucleotide))

def mutated_sequence_based_on_mutation_rate_p(sequence,p_mutation_rate):
    """This function taks in a sequence to be mutated based on the mutation rate p, which is another argument of the function. 
    The function loops through each position of the sequence,
    decides whether the nucleotide at that position is mutated, 
    and adds on to the newly mutated sequence, which will then be return in the end."""
    mutated_sequence = ''
    for nt_position in sequence:
        if mutate_position_based_on_mutation_rate_p(p_mutation_rate):
            mutate_with = random.choice(mutate_with_nucleotides(nt_position))
            mutated_sequence+=mutate_with
        else:
            mutated_sequence+=nt_position
    return(mutated_sequence)

def get_coding_sequence_from_nucleotide_sequence(nt_sequence):
    """This functioin returns the coding sequence as a list when given a nucleotidesequence.
    nt_sequence: nucleoeitde sequence, preferably a protein coding gene sequence"""
    return(list(sliced(nt_seq[frame-1:].upper(),3)))

def translate_coding_sequence(cds_seq):
    """Function translates open reading frames using the coding sequence.
    cds_seq: a list of codons for a given sequemce"""
    codontab=codon_table()
    translation=''
    for codon in cds_seq:
        if codon in codontab:
            translation+=codontab[codon]
    return(translation)


