import pandas as pd
from fmh_omega import dnds
k=15
kdna=k*3
wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/marine_bacteria'
file1=pd.read_csv(f'{wd}/ani.dna.45.extracted.csv').set_index(['A','B'])[['ANI']]
file2=pd.read_csv(f'{wd}/aai.prot.15.extracted.csv').set_index(['A','B'])[['AAI']]
file3=pd.read_csv(f'{wd}/nt_containment15.csv').rename(columns={'containment':'DNA_max_Cfrac'}).set_index(['A','B'])[['DNA_max_Cfrac']]
file4=pd.read_csv(f'{wd}/prot_containment15.csv').rename(columns={'containment':'AA_max_Cfrac'}).set_index(['A','B'])[['AA_max_Cfrac']]
file5=pd.read_csv(f'{wd}/dnds_constant_15.csv').set_index(['A','B']).rename(columns={'dNdS_ratio_constant':'dN/dS'})[['dN/dS']]
data=pd.concat([file1,file2,file3,file4,file5],axis=1)
data['ANI_approx']=data.apply(lambda row: dnds.ANI_approx(row['DNA_max_Cfrac'], kdna), axis=1)
data['AAI_approx']=data.apply(lambda row: dnds.AAI_approx(row['AA_max_Cfrac'], k), axis=1)
data.to_csv(f'{wd}/average_identities_data.k15.csv')