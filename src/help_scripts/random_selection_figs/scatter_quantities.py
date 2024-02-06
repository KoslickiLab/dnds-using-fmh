import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

p = 0.001
select =''
quantity = 'dNdS'
k=7

ground_truth_file = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/ground_truth/{select}_{p}.csv",sep=",",header=0).dropna()
#ground_truth_file = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/real_data_{p}_ksizes_5_30_10002_cdn/ground_truth_{p}.csv",sep=",",header=0)

fmh_file = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds/{select}_selection/dnds_constant_{k}.csv",sep=",",header=0)

if quantity == 'AAI':
    fmh_quantity = 'containment_protein'
elif quantity == 'ANI':
    fmh_quantity = 'containment_nt'
elif quantity == 'PdN':
    fmh_quantity = 'dN'
elif quantity == 'PdS':
    fmh_quantity = 'dS'
elif quantity == 'dNdS':
    fmh_quantity = 'dNdS_ratio_constant'
    ground_truth_file = ground_truth_file[(ground_truth_file['dNdS'] > 0) & (ground_truth_file['dNdS'] < 2)].set_index('Unnamed: 0')[['dNdS']]
    fmh_file = fmh_file[(fmh_file['dNdS_ratio_constant'] > 0) & (fmh_file['dNdS_ratio_constant'] < 2)].set_index('B')[['dNdS_ratio_constant']]
    fmh_file = fmh_file.set_index('B')[['dNdS_ratio_constant']]
    df = pd.concat([ground_truth_file,fmh_file],axis=1).dropna()

#create figure    
plt.scatter(df[quantity], df[fmh_quantity])
plt.title(f'{quantity} range 0-2 {select} {p} 10,0002nt n={len(df)}')
plt.xlabel(f'ground truth')
plt.ylabel(f'fmh k_aa={k}')
tmp = [min(df[quantity]),max(df[quantity])]
plt.plot(tmp,tmp,linestyle="--")
plt.savefig(f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/ground_truth/compared_to_FMH_ground_truth_{quantity}_{select}_{k}_{p}_range0-2.png')
