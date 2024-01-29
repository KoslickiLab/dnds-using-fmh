import pickle
import pandas as pd

""" deprecated function
def grab_containment_from_mat(mat_df,ksize):
    #This function converts matrix into df removes pairwise information
    #When running sourmash comapre, a matrix via a csv file called compare.csv is produced
    mat = pd.read_csv(mat_df,sep=',')
    df = pd.DataFrame(columns=['A','B','containment','ksize'],dtype=int)
    for i in range(len(mat)):
        for j in range(len(mat)):
            if i > j:
                containment = mat.iloc[j,i]
                df = pd.concat([df, pd.DataFrame.from_records([{'A':mat.columns[i],
                'B':mat.columns[j],
                'containment':containment,
                'ksize':ksize}])],join='inner')
    return(df)
""" 

def grab_containment_from_mat_ground_truth(mat_df,ksize):
    """This function converts matrix into df removes pairwise information"""
    """When running sourmash comapre, a matrix via a csv file called ccdompare.csv is produced"""

    #read in df
    #total_sigs=9
    df = pd.read_csv(mat_df,sep=',')

    #record gene header into list
    gene_name_header_list = df.T.index.to_list()

    subset_number = int(len(gene_name_header_list)/2)
    subset = df.iloc[0:subset_number, 0:subset_number]
    
    #make the gene header list into a column to set as index
    subset['A'] = gene_name_header_list[:subset_number]
    
    #create df
    subset = subset.set_index('A').stack().reset_index().rename(columns={'level_1':'B',0:'containment'})
    
    #dont forget to add ksize column!
    subset['ksize']=ksize
    
    return(subset)

def label_selection_pressure(row):
   if row['dN/dS'] == 1:
      return 'neutral'
   if row['dN/dS'] > 1:
      return 'positive'
   if row['dN/dS'] < 1:
      return 'negative'

def change_column_names(csv_file):
    df = pd.read_csv(csv_file,sep=',')
    genome_set = set(df['A'].tolist() + df['B'].tolist())
    for i in genome_set:
        for j in genome_set:
            if i != j:
                df.loc[(df['A'] == j) & (df['B'] == i), ['A', 'B']] = [i, j]
    return(df)

""" deprecated function
def grab_max_nt_containment_from_containment_csv_file(csv_file):
    #Return the max C(A,B)
    data = pd.read_csv(csv_file,sep=',')
    genome_set = set(data['A'].tolist() + data['B'].tolist())
    data = data.set_index(['A','B'])
    for i in genome_set:
        for j in genome_set:
            if (i,j) in data.index.tolist() and (j,i) in data.index.tolist():
                if data.loc[i].loc[j]['containment_nt'] > data.loc[j].loc[i]['containment_nt']:
                    data = data.drop(index=(j, i))
                else:
                    data =data.drop(index=(i,j))
    data['selection_pressure'] = data.apply(label_selection_pressure, axis=1)
    return(data)
"""

""" deprecated function
def grab_max_protein_containment_from_containment_csv_file(csv_file):
    #Return the max C(A,B)
    data = pd.read_csv(csv_file,sep=',')
    genome_set = set(data['A'].tolist() + data['B'].tolist())
    data = data.set_index(['A','B'])
    for i in genome_set:
        for j in genome_set:
            if (i,j) in data.index.tolist() and (j,i) in data.index.tolist():
                if data.loc[i].loc[j]['containment_protein'] > data.loc[j].loc[i]['containment_protein']:
                    data = data.drop(index=(j, i))
                else:
                    data =data.drop(index=(i,j))
    data['selection_pressure'] = data.apply(label_selection_pressure, axis=1)
    return(data)
"""