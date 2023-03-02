"""Module containing code for calculating dNdS from min hash containment indexes"""

import pandas as pd

def calc_Pnt(nt_containment,k):
    """ 
    Returns nucleotide mutation rate between two protein-coding sequences (or any other type of nucleotide sequences)
    nt_containment: The containment index between two nucleotide sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    nt_k = 3*k
    return(1 - nt_containment**(1/nt_k))

def calc_PdN(protein_containment,k):
    """
    Returns nonsynonymous mutation rate between two protein sequences
    protein_containment: The containment index between two protein sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    return(1-protein_containment**(1/k))

def calc_nomutation(nt_containment,k):
    """
    Returns no mutation rate between two protein-coding sequences (or any other type of nucleotide sequences).
    Uses calc_Pnt() function for estimation.
    nt_containment: The containment index between two nucleotide sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    return((1-calc_Pnt(nt_containment,k))**3)

def calc_PdS(protein_containment,nt_containment,k):
    """
    Returns synonymous mutation rate between two protein sequences.
    Uses calc_nomutation() function for estimation.
    nt_containment: The containment index between two nucleotide sequences (this is a float)
    protein_containment: The containment index between two protein sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    PdS = 1-calc_PdN(protein_containment,k) - (calc_nomutation(nt_containment,k))
    return(PdS)

def dNdS_ratio(protein_containment,nt_containment,k):
    """
    Returns dN/dS ratio between two protein sequences.
    Uses calc_PdN() and calc_PdS() functions for estimation
    nt_containment: The containment index between two nucleotide sequences (this is a float)
    protein_containment: The containment index between two protein sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    return((calc_PdN(protein_containment,k))/(calc_PdS(protein_containment,nt_containment,k)))

def report_dNdS(nt_containment_df,prot_containment_df):
    """
    Returns dataframe of dNdS reports between all pairwise estimations of protein coding sequences.
    Uses dNdS_ratio() function to estimate dN/dS ratio.
    nt_containment_df: dataframe with prepared containment indexes of nucleotide sequences that is 
        produced by grab_containment_from_mat() function from reportCI.py. 
        Dataframe contains ref, query, containment index, and ksize.
    prot_containment_df: dataframe with prepared containment indexes of protein sequences that is 
        produced by grab_containment_from_mat() function from reportCI.py. 
        Dataframe contains ref, query, containment index, and ksize.
    """
    #read in nt_containment dataframe file and change column names
    nt_df = nt_containment_df.rename(columns={'containment':'containment_nt'})    
    #read in nt_containment dataframe file and change column names
    protein_df = prot_containment_df.rename(columns={'containment':'containment_protein'})
    #join df into one
    df = pd.merge(nt_df, protein_df, on=['A','B','ksize'])
    #apply function
    df['dNdS_ratio'] = dNdS_ratio(nt_containment=df['containment_nt'],protein_containment=df['containment_protein'],k=df['ksize'])

    #report
    return(df)
