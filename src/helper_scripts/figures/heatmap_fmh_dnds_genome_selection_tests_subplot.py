import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#create figure that represents median Cfracs across ksizes for both DNA and Protein
clmns=['dN/dS']

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k7_10scale_cores1000/fmh_omega_7.csv'
file1=pd.read_csv(f'{file}',sep=',')[clmns]
print(file1.describe())

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k7_10scale_cores1000/fmh_omega_7.csv'
file2=pd.read_csv(f'{file}',sep=',')[clmns]
print(file2.describe())

# Create subplots
fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 1]})
max_value = max(max(file1['dN/dS']),max(file2['dN/dS']))

# Plot the first heatmap
heatmap1 = sns.heatmap(file1, ax=axes[0], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap='coolwarm', vmin=0, vmax=max_value, center=1,
                cbar_kws={'label': 'FMH dN/dS'},cbar=False)
axes[0].set_title('Positive')
axes[0].set_xticks([])

heatmap2 = sns.heatmap(file2, ax=axes[1], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap='coolwarm', vmin=0, vmax=max_value, center=1,
                cbar_kws={'label': 'FMH dN/dS'}, cbar=False)
axes[1].set_title('Negative')
axes[1].set_xticks([])

# Create a common color bar
cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])  # Adjust the position and size as needed
cbar = fig.colorbar(heatmap2.get_children()[0], cax=cbar_ax, orientation='vertical', label='FMH dN/dS')


# Adjust layout
#plt.tight_layout()
plt.subplots_adjust(wspace=0.03)
#fig.text(0.5, 0.03, 'FMH dN/dS', ha='center', va='center')

plt.suptitle('Selection simulation on E. coli protein rep genome, p=0.01,ksize=7,scale=10')
#fig.figure.savefig(f"{wd}/ANI_and_AAI.png",bbox_inches='tight')
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/fmh_dnds_genome_selection_tests_heatmaps.tiff",bbox_inches='tight') 
