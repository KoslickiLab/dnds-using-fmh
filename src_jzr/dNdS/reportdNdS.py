"""Module containing code for calculating dNdS from min hash containment indexes"""

import pandas as pd

def calc_Pnt(nt_containment,k):
    return(1 - (nt_containment)**(1/k))

def calc_PdN(protein_containment,k):
    return(1-(protein_containment)**(1/k))

def calc_nomutation(nt_containment,k):
    return((1-calc_Pnt(nt_containment,k))**3)

def calc_PdS(protein_containment,nt_containment,k):
    return(1 - calc_nomutation(nt_containment,k))

def dNdS_ratio(protein_containment,nt_containment,k):
    return(calc_PdN(protein_containment,k)/calc_PdS(protein_containment,nt_containment,k))

def report_dNdS(nt_containment_df,prot_containment_df):
    #read in nt_containment csv file and change column names
    nt_df = nt_containment_df.rename(columns={'containment':'containment_nt'})
    
    #read in nt_containment csv file and change column names
    protein_df = prot_containment_df.rename(columns={'containment':'containment_protein'})

    #join df into one
    df = pd.merge(nt_df, protein_df, on=['A','B','ksize'])

    #apply function
    df['dNdS_ratio'] = dNdS_ratio(df['containment_nt'],df['containment_protein'],df['ksize'])

    #report
    return(df)
