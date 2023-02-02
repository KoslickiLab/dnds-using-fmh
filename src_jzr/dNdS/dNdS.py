"""Module containing code for calculating dNdS from min hash containment indexes"""

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

