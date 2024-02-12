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

