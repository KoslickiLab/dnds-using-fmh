import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr
from matplotlib.ticker import FormatStrFormatter


p = 0.01
select ='positive'
len=10002
d='dS'

if d == 'dN':
    col='Ka'
elif d == 'dS':
    col='Ks'

file1 = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/ground_truth/ground_truth_0.01_positive.csv",sep=",",header=0)
d_file1 = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/positive_selection_queries_10002_0.01.axt.kaks",sep="\t",header=0)

file2 = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/ground_truth/ground_truth_0.01_negative.csv",sep=",",header=0)
d_file2 = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/positive_selection_queries_10002_0.01.axt.kaks",sep="\t",header=0)

# Create subplots with adjusted layout
fig, axes = plt.subplots(2, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 1]})

#create figure    
axes[0, 0].scatter(file1["ANI"],d_file1[col], color='red')
axes[0, 0].set_title('Positive')
axes[0, 0].set_ylabel(f'{d}')


axes[0, 1].scatter(file1["AAI"],d_file1[col], color='red')
axes[0, 1].set_title('Positive')

axes[1, 0].scatter(file2["ANI"],d_file2[col], color='blue')
axes[1, 0].set_title('Negative')
axes[1, 0].set_ylabel(f'{d}')
axes[1, 0].set_xlabel('ANI')
axes[1, 0].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))


axes[1, 1].scatter(file2["AAI"],d_file2[col], color='blue')
axes[1, 1].set_title('Negative')
axes[1, 1].set_xlabel('AAI')

axes[0, 0].text(-0.08, 1.1, f'({chr(ord("a") + 0)})', transform=axes[0, 0].transAxes,
                 size=12, weight='bold')
axes[1, 0].text(-0.08, 1.1, f'({chr(ord("a") + 1)})', transform=axes[0, 1].transAxes,
                 size=12, weight='bold')
axes[0, 1].text(-0.08, 1.1, f'({chr(ord("a") + 2)})', transform=axes[1, 0].transAxes,
                 size=12, weight='bold')
axes[1, 1].text(-0.08, 1.1, f'({chr(ord("a") + 3)})', transform=axes[1, 1].transAxes,
                 size=12, weight='bold')

plt.subplots_adjust(hspace=0.4)

plt.savefig(f'/data/jzr5814/sourmash_dnds_estimation/thesis_figures/scatter_ground_truth_AXI_to_{d}_{p}.png')
