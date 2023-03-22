import pandas as pd
import subprocess

file='/data/jzr5814/repositories/dnds-using-fmh/src/sourmash_dnds_estimation/merged_KO_gene_seq_lens.txt'
data=pd.read_csv(file,delimiter='\t',skiprows=[0],header=None).rename(columns={0:"gene",1:"nt_length"}).sort_values(by="nt_length")    

fna_output = open('real_protein_coding_genes.fna','w')

for i in [[500,600],[1000,2000],[2000,3000],[4000,5000],[8000,9000],[10000,11000]]:
    gene=data[(i[0] <= data['nt_length']) & (data['nt_length'] < i[1])].iloc[0,0]
    sequence=data[(i[0] <= data['nt_length']) & (data['nt_length'] < i[1])].iloc[0,2]
    fna_output.write(f'>{gene}_{len(sequence)}\n')
    fna_output.write(f'{sequence}\n')


""""gene_500nt=data[(500 <= data['nt_length']) & (data['nt_length'] < 600)].iloc[0,0]
sequence_500nt=data[(500 <= data['nt_length']) & (data['nt_length'] < 600)].iloc[0,2]

gene_1000nt=data[(1000 <= data['nt_length']) & (data['nt_length'] < 2000)].iloc[0,0]
sequence_1000nt=data[(1000 <= data['nt_length']) & (data['nt_length'] < 2000)].iloc[0,2]

gene_2000nt=data[(2000 <= data['nt_length']) & (data['nt_length'] < 3000)].iloc[0,0]
sequence_2000nt=data[(2000 <= data['nt_length']) & (data['nt_length'] < 3000)].iloc[0,2]

gene_4000nt=data[(4000 <= data['nt_length']) & (data['nt_length'] < 5000)].iloc[0,0]
sequence_4000nt=data[(4000 <= data['nt_length']) & (data['nt_length'] < 5000)].iloc[0,2]

gene_8000nt=data[(8000 <= data['nt_length']) & (data['nt_length'] < 9000)].iloc[0,0]
sequence_8000nt=data[(8000 <= data['nt_length']) & (data['nt_length'] < 9000)].iloc[0,2]

gene_10000nt=data[(10000 <= data['nt_length']) & (data['nt_length'] < 11000)].iloc[0,0]
sequence_10000nt=data[(10000 <= data['nt_length']) & (data['nt_length'] < 11000)].iloc[0,2]
"""
