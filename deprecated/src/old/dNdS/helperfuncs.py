import pickle
import pandas as pd
from Bio import SeqIO

def containments(mat_df,ksize):
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
    #create df
    subset = subset.set_index('A').stack().reset_index().rename(columns={'level_1':'B',0:'containment'})
    #dont forget to add ksize column!
    subset['ksize']=ksize
    return(subset)

def translate_CDS(cds_fasta, output_filename):
    """User has input fasta file with CDS sequences of a genome"""
    """This function translates each CDS found in the FASTA file"""
    sequences = SeqIO.parse(open(cds_fasta),'fasta')
    with open(output_filename) as out_file:
        for cds in sequences:
            name, cds_seq = cds.id, str(cds.seq)
            aa_seq = Seq(cds_seq).translate()
            out_file.write(''.join['>',name,'\n'])
            out_file.write(''.join[aa_seq,'\n'])
        out_file.close()
