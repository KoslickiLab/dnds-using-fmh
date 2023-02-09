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
