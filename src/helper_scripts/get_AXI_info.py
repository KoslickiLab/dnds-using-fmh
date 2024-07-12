from fmh_omega import helperfuncs, dnds
import pandas as pd

wd='/data/jzr5814/sourmash_dnds_estimation/tests/results/genomic_dnds/marine_bacteria'

ani_mat_df = helperfuncs.extract_axi_matrix(f'{wd}/ani.dna.21.csv','ANI').set_index(['A','B'])

aai_mat_df = helperfuncs.extract_axi_matrix(f'{wd}/aai.prot.7.csv',"AAI").set_index(['A','B'])

dnds_df = pd.read_csv(f'{wd}/dnds_constant_7.csv').set_index(['A','B']).rename(columns={"dNdS_ratio_constant":"dN/dS",'containment_protein':'AA_max_Cfrac','containment_nt':'DNA_max_Cfrac'})[['DNA_max_Cfrac','AA_max_Cfrac','dN/dS']]
dnds_df['ANI_approx'] = dnds_df.apply(lambda row: dnds.ANI_approx(row['DNA_max_Cfrac'], 21), axis=1)
dnds_df['AAI_approx'] = dnds_df.apply(lambda row: dnds.AAI_approx(row['AA_max_Cfrac'], 7), axis=1)

data = pd.concat([ani_mat_df, aai_mat_df, dnds_df], axis=1, join="inner")
data.to_csv(f'{wd}/average_identities_data.k7.csv')
