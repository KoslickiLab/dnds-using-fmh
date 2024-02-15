import pickle
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import os

def extract_filename_without_extension(file_path):
    #return file_path.split('/')[-1].split('.')[0]
    return file_path.split('/')[-1].split('.')[0]

def containments(mat_df,ksize,multiple):
    """This function converts matrix into df removes pairwise information"""
    """When running sourmash compare, a matrix via a csv file is produced"""
    #read in df
    df = pd.read_csv(mat_df,sep=',')
    #record gene header into list
    gene_name_header_list = df.T.index.to_list()
    subset_number = int(len(gene_name_header_list)/2) #containment matrix produced by sourmash is not perfect square
    subset = df.iloc[0:subset_number, 0:subset_number]
    #make the gene header list into a column to set as index
    subset['A'] = gene_name_header_list[:subset_number]
    subset = subset.set_index('A').stack().reset_index().rename(columns={'level_1':'B',0:'containment'})
    if multiple=='yes':
        subset['A']=subset['A'].apply(extract_filename_without_extension)
        print('change A')
        print(subset['A'])
        subset['B']=subset['B'].apply(extract_filename_without_extension)
        print('change B')
        print(subset['B'])
        
    #dont forget to add ksize column!
    subset['ksize']=ksize
    return(subset)

def divisible_by_3(sequence):
    #makes sure a sequence is divisible by 3 for translation
    if len(sequence) % 3 != 0:
        if len(sequence+'N') % 3 != 0:
            return(sequence+'NN')
        else:
            return(sequence+'N')
    else:
        return(sequence)

def translate_CDS(cds_fasta, out_name):
    """User has input fasta file with CDS sequences of a genome"""
    """This function translates each CDS found in the FASTA file"""
    sequences = SeqIO.parse(open(cds_fasta),'fasta')
    with open(f'{out_name}','w') as out_file:
        for cds in sequences:
            name, nt_seq = cds.id, str(cds.seq)
            cds_seq = divisible_by_3(nt_seq)
            aa_seq = str(Seq(cds_seq).translate())
            out_file.write(''.join(['>',name,'\n']))
            out_file.write(''.join([aa_seq,'\n']))
        out_file.close()

def return_protein_klist_parameters(kmer_list):
    sm_klist = ',k='.join(kmer_list.split(','))
    return(sm_klist)

def return_dna_klist_parameters(kmer_list):
    temp_klist = kmer_list.split(',')
    sm_klist=[]
    for k in temp_klist:
        sm_klist.append(str(int(k)*3))
    return(',k='.join(sm_klist))

def return_signature_list(working_dir, molecule):
    sigs = []
    for file in os.listdir(f"{working_dir}/signatures"):
        if file.endswith(f".{molecule}.sig.gzip"):
            sigs.append(file)
    return sigs

#The following functions create a random sequence for selection

def get_coding_sequence_from_nucleotide_sequence(nt_sequence):
    """This functioin returns the coding sequence as a list when given a nucleotide sequence.
    nt_sequence: nucleotide sequence, preferably a protein coding gene sequence"""
    coding_sequence = list(sliced(nt_sequence,3))     #coding sequence
    codons_to_remove = ['TGA', 'TAA', 'TAG']     #list of stop codons
    filtered_codons = [codon for codon in coding_sequence if codon not in codons_to_remove]     # Remove stop codons
    return(filtered_codons)

def mutate_position_based_on_mutation_rate_p(p_mutation_rate):
    """This function returns a boolean probability of a mutation based off a given mutation rate p
    p_mutation_rate = 1 - Cfrac(a,B)**(1/k)
    Cfrac is the containment index between two sequences"""
    p = random.random()
    if p <= p_mutation_rate:
        return(True)
    else:
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

def mutate_position(nt_position, p_mutation_rate):
    prate = mutate_position_based_on_mutation_rate_p(p_mutation_rate)
    if prate:
        mutate_with = random.choice(mutate_with_nucleotides(nt_position.upper()))
        return(mutate_with)
    else:
        return(nt_position)

def positive_selection_outcome(codon,p_mutation_rate):
    mutated_codon = ''
    if random.choice([0, 1, 2])==random.choice([0, 1, 2]): 
        for position in range(len(codon)):
            nt_position = codon[position]
            if position == 0 or position == 1:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
            elif position == 2:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
    else:
        for position in range(len(codon)):
            nt_position = codon[position]
            if position == 0 or position == 1:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
            else:
                mutated_codon+=nt_position
    return(mutated_codon)

def negative_selection_outcome(codon,p_mutation_rate):
    mutated_codon = ''
    if random.choice([0, 1, 2])==random.choice([0, 1, 2]): 
        for position in range(len(codon)):
            nt_position = codon[position]
            if position == 0 or position == 1:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
            elif position == 2:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
    else:
        for position in range(len(codon)):
            nt_position = codon[position]
            if position == 2:
                mutated_codon += mutate_position(nt_position, p_mutation_rate)
            else:
                mutated_codon+=nt_position
    return(mutated_codon)

def positive_selection_based_on_mutation_rate_p(sequence,p_mutation_rate):
    """This function takes in a sequence to be mutated based on the mutation rate p, which is another argument of the function. 
    The function loops through each position of the sequence,
    decides whether the nucleotide at that position is mutated, 
    and adds on to the newly mutated sequence, which will then be return in the end."""
    mutated_sequence = ''
    coding_sequence = list(sliced(sequence,3))
    for codon in coding_sequence:
        mutated_sequence+=positive_selection_outcome(codon.upper(), p_mutation_rate)
    return(mutated_sequence.upper())

def negative_selection_based_on_mutation_rate_p(sequence,p_mutation_rate):
    """This function takes in a sequence to be mutated based on the mutation rate p, which is another argument of the function. 
    The function loops through each position of the sequence,
    decides whether the nucleotide at that position is mutated, 
    and adds on to the newly mutated sequence, which will then be return in the end."""
    mutated_sequence = ''
    coding_sequence = list(sliced(sequence,3))
    for codon in coding_sequence:
        mutated_sequence+=negative_selection_outcome(codon.upper(), p_mutation_rate)
    return(mutated_sequence.upper())
