import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics

#input options

p1=0.1
p2=0.01
p3=0.001

wd1 = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p1}/kaks_NG'
wd2 = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p2}/kaks_NG'
wd3 = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p3}/kaks_NG'

data_positive1 = f'positive_selection_queries_10002_{p1}.axt.kaks'
data_positive2 = f'positive_selection_queries_10002_{p2}.axt.kaks'
data_positive3 = f'positive_selection_queries_10002_{p3}.axt.kaks'

positive1 = pd.read_csv(f'{wd1}/{data_positive1}',sep='\t')
positive2 = pd.read_csv(f'{wd2}/{data_positive2}',sep='\t')
positive3 = pd.read_csv(f'{wd3}/{data_positive3}',sep='\t')

positive1 = positive1['Ka/Ks'].dropna().tolist()
positive2 = positive2['Ka/Ks'].dropna().tolist()
positive3 = positive3['Ka/Ks'].dropna().tolist()

data_negative1= f'negative_selection_queries_10002_{p1}.axt.kaks'
data_negative2= f'negative_selection_queries_10002_{p2}.axt.kaks'
data_negative3= f'negative_selection_queries_10002_{p3}.axt.kaks'

negative1 = pd.read_csv(f'{wd1}/{data_negative1}',sep='\t')
negative2 = pd.read_csv(f'{wd2}/{data_negative2}',sep='\t')
negative3 = pd.read_csv(f'{wd3}/{data_negative3}',sep='\t')

negative1 = negative1['Ka/Ks'].dropna().tolist()
negative2 = negative2['Ka/Ks'].dropna().tolist()
negative3 = negative3['Ka/Ks'].dropna().tolist()

# Set the Seaborn style and remove the top and right spines
sns.set_style("whitegrid")

# Create subplots
fig, axes = plt.subplots(3, 1, figsize=(12, 5), sharey=True, sharex=True)

# plotting
sns.histplot(positive1, label='Positive', color='red', linewidth=0, ax=axes[0])
sns.histplot(negative1, label='Negative', color='blue', linewidth=0, ax=axes[0])
#sns.despine(ax=axes[0], top=True, right=True)
axes[0].set_title(f'{p1}', fontsize=15)
axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axes[0].grid(axis='y', linestyle='--')

# plotting
sns.histplot(positive2, label='Positive', color='red', linewidth=0, ax=axes[1])
sns.histplot(negative2, label='Negative', color='blue', linewidth=0, ax=axes[1])
#sns.despine(ax=axes[1], top=True, right=True)
axes[1].set_title(f'{p2}', fontsize=15)
axes[1].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axes[1].grid(axis='y', linestyle='--')

# plotting
sns.histplot(positive3, label='Positive', color='red', linewidth=0, ax=axes[2])
sns.histplot(negative3, label='Negative', color='blue', linewidth=0, ax=axes[2])
#sns.despine(ax=axes[2], right=True)
axes[2].set_title(f'{p3}', fontsize=15)
axes[2].set_xlabel("NG86-Based dN/dS")

# Add legend to the first plot
axes[0].legend()


axes[0].text(-0.08, 1.1, f'({chr(ord("a") + 0)})', transform=axes[0].transAxes,
                 size=12, weight='bold')
axes[1].text(-0.08, 1.1, f'({chr(ord("a") + 1)})', transform=axes[1].transAxes,
                 size=12, weight='bold')
axes[2].text(-0.08, 1.1, f'({chr(ord("a") + 2)})', transform=axes[2].transAxes,
                 size=12, weight='bold')

# Remove grid from both subplots
for ax in axes:
    ax.grid(False)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 35)
    ax.axvline(x=1, color='grey', linestyle='--', label='Vertical Line at x=3')

axes[2].set_xticks([0, 1, 2, 4, 6, 8, 10, 12])

#plot_title = f"NG86-Based dN/dS of Selectively Mutated Sequences of Length 10,002 nt"
#plt.suptitle(plot_title,x=0.5,y=1)

# Set the Seaborn style and remove the top and right spines

plt.subplots_adjust(hspace=0.4)

plt.savefig(f'/data/jzr5814/sourmash_dnds_estimation/thesis_figures/kaks_selection_subplot_dist_thesis.png',bbox_inches='tight')
