import pickle
import pandas as pd

def grab_containment_from_mat(mat_df,ksize):
    """This function converts matrix into df removes pairwise information"""
    """When running sourmash comapre, a matrix via a csv file called compare.csv is produced"""
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

def grab_containment_from_mat_ground_truth(mat_df,ksize):
    """This function converts matrix into df removes pairwise information"""
    """When running sourmash comapre, a matrix via a csv file called compare.csv is produced"""
    #df = pd.read_csv(mat_df,sep=',').T.reset_index()[['index',0]].rename(columns={'index':'B',0:'containment'})
    #df['A']='ref_gene'
    #df['ksize']=ksize

    #read in df
    df = pd.read_csv(mat_df,sep=',')

    #record gene header into list
    gene_name_header_list = df.T.index.to_list()

    #make the gene header list into a column to set as index
    df['A'] = gene_name_header_list
    
    #create df
    df = df.set_index('A').stack().reset_index().rename(columns={'level_1':'B',0:'containment'})
    
    #dont forget to add ksize column!
    df['ksize']=ksize

    return(df)
