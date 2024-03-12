"""Module containing code for calculating dNdS from min hash containment indexes"""

import pandas as pd
import logging
from loguru import logger
import time

def calc_PdN(protein_containment,k):
    """
    Returns nonsynonymous mutation rate between two protein sequences
    protein_containment: The containment index between two protein sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    return(1-protein_containment**(1/k))

def calc_PdS(protein_containment,nt_containment,k):
    """
    Returns synonymous mutation rate between two protein sequences.
    Uses calc_nomutation() function for estimation.
    nt_containment: The containment index between two nucleotide sequences (this is a float)
    protein_containment: The containment index between two protein sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    nt_k = 3*k
    P_nt_mut = 1 - nt_containment**(1/nt_k)
    P_no_mut = (1-P_nt_mut)**3
    PdS = 1 - calc_PdN(protein_containment,k) - P_no_mut
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

def dNdS_ratio_with_constant(protein_containment,nt_containment,k):
    """
    Returns dN/dS ratio between two protein sequences.
    Uses calc_PdN() and calc_PdS() functions for estimation
    nt_containment: The containment index between two nucleotide sequences (this is a float)
    protein_containment: The containment index between two protein sequences (this is a float)
    k: Identify the ksize used to produce containment index (this is an integer)
    """
    constant = 0.77/2.23
    try:
        logger.info(f"Estimating dN/dS using FMH")
        dNdS = calc_PdN(protein_containment,k)/calc_PdS(protein_containment,nt_containment,k)
        dNdS_constant = dNdS*constant
        logger.success(f"Successfully estimated")
    except ZeroDivisionError as e:
        logging.error("ZeroDivisionError: ignore undefined dN/dS estimation")
    return(dNdS_constant)

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
    #read in nt_containment and protein_containment dataframe file and change column names
    nt_df = nt_containment_df.rename(columns={'containment':'DNA_Cfrac'})
    protein_df = prot_containment_df.rename(columns={'containment':'AA_Cfrac'})
    #join df into one
    df = pd.merge(nt_df, protein_df, on=['A','B','ksize'])
    #apply function
    df['PdN'] = (calc_PdN(protein_containment=df['AA_Cfrac'],k=df['ksize']))
    df['PdS'] = (calc_PdS(protein_containment=df['AA_Cfrac'],nt_containment=df['DNA_Cfrac'],k=df['ksize']))
    df['PdN/PdS'] = dNdS_ratio(nt_containment=df['DNA_Cfrac'],protein_containment=df['AA_Cfrac'],k=df['ksize'])
    df['dN/dS'] = dNdS_ratio_with_constant(nt_containment=df['DNA_Cfrac'],protein_containment=df['AA_Cfrac'],k=df['ksize'])

    #report
    return(df)

def report_dNdS_multisearch(dna_cfrac_csv,protein_cfrac_csv,ksize):
    """
    Returns dataframe of dNdS reports between all pairwise estimations of protein coding sequences.
    Uses dNdS_ratio() function to estimate dN/dS ratio.
    dna_cfrac_csv: multisearch results csv file for dna 
        Header of csv file contains ref, query, containment index, and ksize.
    protein_cfrac_csv: multisearch results csv file for protein
        Header of csv file contains ref, query, containment index, and ksize.
    ksize: ksize of protein
    """
    #read in nt_containment and protein_containment dataframe file and change column names
    #nt_df = nt_containment_df.rename(columns={'containment':'DNA_Cfrac'})
    dna_cfrac = pd.read_csv(f'{dna_cfrac_csv}',sep=",")[['query_name','match_name','containment']].rename(columns={'query_name':'A','match_name':'B','containment':'DNA_Cfrac'})
    #protein_df = prot_containment_df.rename(columns={'containment':'AA_Cfrac'})
    protein_cfrac = pd.read_csv(f'{protein_cfrac_csv}',sep=",")[['query_name','match_name','containment']].rename(columns={'query_name':'A','match_name':'B','containment':'AA_Cfrac'})

    #join df into one
    #df = pd.merge(nt_df, protein_df, on=['A','B','ksize'])
    merge_df = pd.merge(dna_cfrac, protein_cfrac, on=['A','B'])
    merge_df['ksize'] = int(ksize)

    #apply function
    merge_df['PdN'] = (calc_PdN(protein_containment=merge_df['AA_Cfrac'],k=merge_df['ksize']))
    merge_df['PdS'] = (calc_PdS(protein_containment=merge_df['AA_Cfrac'],nt_containment=merge_df['DNA_Cfrac'],k=merge_df['ksize']))
    merge_df['PdN/PdS'] = dNdS_ratio(nt_containment=merge_df['DNA_Cfrac'],protein_containment=merge_df['AA_Cfrac'],k=merge_df['ksize'])
    merge_df['dN/dS'] = dNdS_ratio_with_constant(nt_containment=merge_df['DNA_Cfrac'],protein_containment=merge_df['AA_Cfrac'],k=merge_df['ksize'])

    #report
    return(merge_df)


def report_dNdS_pairwise(dna_cfrac_csv,protein_cfrac_csv,ksize):
    """
    Returns dataframe of dNdS reports between all pairwise estimations of protein coding sequences.
    it seems that the max containment is used when reporting pairwise containments 
    Uses dNdS_ratio() function to estimate dN/dS ratio.
    dna_cfrac_csv: multisearch results csv file for dna 
        Header of csv file contains ref, query, containment index, and ksize.
    protein_cfrac_csv: multisearch results csv file for protein
        Header of csv file contains ref, query, containment index, and ksize.
    ksize: ksize of protein
    """
    #read in nt_containment and protein_containment dataframe file and change column names
    #nt_df = nt_containment_df.rename(columns={'containment':'DNA_Cfrac'})
    start_time = time.time()
    dna_cfrac = pd.read_csv(f'{dna_cfrac_csv}',sep=",")[['query_name','match_name','max_containment']].rename(columns={'max_containment':'DNA_max_Cfrac'})
    dna_cfrac['A,B'] = dna_cfrac[['query_name', 'match_name']].apply(sorted, axis=1).apply(tuple)
    #protein_df = prot_containment_df.rename(columns={'containment':'AA_Cfrac'})
    protein_cfrac = pd.read_csv(f'{protein_cfrac_csv}',sep=",")[['query_name','match_name','max_containment']].rename(columns={'max_containment':'AA_max_Cfrac'})
    protein_cfrac['A,B'] = protein_cfrac[['query_name', 'match_name']].apply(sorted, axis=1).apply(tuple)
    #join df into one
    merge_df = pd.merge(dna_cfrac, protein_cfrac, on='A,B', how='inner')
    merge_df['ksize'] = int(ksize)

    #apply function
    merge_df['PdN'] = (calc_PdN(protein_containment=merge_df['AA_max_Cfrac'],k=merge_df['ksize']))
    merge_df['PdS'] = (calc_PdS(protein_containment=merge_df['AA_max_Cfrac'],nt_containment=merge_df['DNA_max_Cfrac'],k=merge_df['ksize']))
    merge_df['PdN/PdS'] = dNdS_ratio(nt_containment=merge_df['DNA_max_Cfrac'],protein_containment=merge_df['AA_max_Cfrac'],k=merge_df['ksize'])
    merge_df['dN/dS'] = dNdS_ratio_with_constant(nt_containment=merge_df['DNA_max_Cfrac'],protein_containment=merge_df['AA_max_Cfrac'],k=merge_df['ksize'])

    end_time = time.time()
    time_logged = end_time-start_time
    print(f"Successfully estimated dN/dS in {time_logged} secconds")
    #report
    return(merge_df)

