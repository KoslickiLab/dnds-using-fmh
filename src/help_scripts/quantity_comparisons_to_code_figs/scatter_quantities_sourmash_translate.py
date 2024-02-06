import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

p = 0.001
quantity = 'dNdS_ratio_constant'
k=7

columns = ['B','A','containment_nt','containment_protein','dN','dS','dNdS_ratio_constant']
ground_truth_file = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p}_ksizes_5_30_10002_cdn/dnds_constant_{k}.csv",sep=",",header=0)[columns]
fmh_file = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p}_ksizes_5_7_10002_cdn_sourmash_translate/dnds_constant_{k}.csv",sep=",",header=0)[columns]
fmh_file = fmh_file[(fmh_file['A'] == 'gene_10002') & (fmh_file['B'] != 'gene_10002')]

if quantity == 'containment_protein':
    prev_quantity = 'prev_AAI'
    trans_quantity = 'trans_AAI'
    title_quantity='AAI'
elif quantity == 'containment_nt':
    prev_quantity = 'prev_ANI'
    trans_quantity = 'trans_ANI'
    title_quantity='ANI'
elif quantity == 'dN':
    prev_quantity = 'prev_PdN'
    trans_quantity = 'trans_PdN'
    title_quantity = quantity
elif quantity == 'dS':
    prev_quantity = 'prev_PdS'
    trans_quantity = 'trans_PdS'
    title_quantity = quantity
elif quantity == 'dNdS_ratio_constant':
    prev_quantity = 'prev_dNdS'
    trans_quantity = 'trans_dNdS'
    ground_truth_file = ground_truth_file[(ground_truth_file[quantity] > 0) & (ground_truth_file[quantity] < 2)].set_index('B')[[quantity]].rename(columns={quantity:prev_quantity})
    #fmh_file = fmh_file[(fmh_file[quantity] > 0) & (fmh_file[quantity] < 2)].set_index('B')[[quantity]].rename(columns={quantity:trans_quantity})
    fmh_file = fmh_file.set_index('B')[[quantity]].rename(columns={quantity:trans_quantity})
    fmh_file = fmh_file[[trans_quantity]]
    quantity='dNdS'
    df = pd.concat([ground_truth_file,fmh_file],axis=1).dropna()

#create figure    
plt.scatter(df[prev_quantity], df[trans_quantity])
#plt.scatter(ground_truth_file[quantity], fmh_file[quantity])
plt.title(f'{quantity} {p} 10,0002nt n={len(df)}')
#plt.title(f'{title_quantity} {p} 10,0002nt')
plt.xlabel(f'previously estimated k_aa={k}')
plt.ylabel(f'sourmash translate estimated k_aa={k}')
tmp = [min(df[prev_quantity]),max(df[prev_quantity])]
#tmp = [min(ground_truth_file[quantity]),max(ground_truth_file[quantity])]
plt.plot(tmp,tmp,linestyle="--")
plt.savefig(f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_0.001_ksizes_5_7_10002_cdn_sourmash_translate/previous_and_sourmash_translate_estimations_{quantity}_{k}_{p}.png')
