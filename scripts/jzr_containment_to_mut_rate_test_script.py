from dnds import *
from helperfuncs import *
import string
import random
from matplotlib import pyplot as plt
import pandas as pd

if __name__ == "__main__":
    k = 7
    correct_dnds_values = []
    estimated_dnds_values_fmh = []

    nt_df=pd.read_csv('/data/jzr5814/repositories/dnds-using-fmh/src/sourmash_dnds_estimation/nt_containment.csv')
    nt_containment_list = nt_df['containment'].tolist()

    prot_df=pd.read_csv('/data/jzr5814/repositories/dnds-using-fmh/src/sourmash_dnds_estimation/prot_containment.csv')
    prot_containment_list = prot_df['containment'].tolist()

    print('approx_dnds_using_fmh')
    for containment in range(len(nt_containment_list)):
        nt_k=3*k
        aa_k=k
        nt_containment=nt_containment_list[containment]
        aa_containemnt=prot_containment_list[containment]
        p_nt = containment_to_mut_rate(nt_containment, nt_k)
        p_aa = containment_to_mut_rate(aa_containemnt, aa_k)

        try:
            dnds = p_aa / ( 1.0 - p_aa - (1.0 - p_nt)**3 )
            print(dnds)
        except:
            if p_nt == p_aa:
                if p_nt == 0.0:
                    #print('No mutations! Identical.')
                    print(-1)
                else:
                    dnds = float('Infinity')
            else:
                #print('Some unknown error occurred!')
                print(-2)
        

