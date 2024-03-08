import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#create figure that represents median Cfracs across ksizes for both DNA and Protein
clmns=['dN/dS']

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k7_1scale_cores1000/fmh_omega_7.csv'
file1=pd.read_csv(f'{file}',sep=',')[clmns]

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k7_10scale_cores1000/fmh_omega_7.csv'
file2=pd.read_csv(f'{file}',sep=',')[clmns]

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/negative_k7_100scale_cores1000/fmh_omega_7.csv'
file3=pd.read_csv(f'{file}',sep=',')[clmns]

# Create subplots
fig, axes = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1, 1, 1]})

# Get the colormap and reverse it
cmap = plt.get_cmap('Blues')
reversed_cmap = cmap.reversed()

# Plot the first heatmap
heatmap1 = sns.heatmap(file1, ax=axes[0], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap=reversed_cmap, vmin=0, vmax=1,
                cbar_kws={'label': 'FMH dN/dS'},cbar=False)
axes[0].set_title('1')
axes[0].set_xlabel('')
print(file1.mean())

heatmap2 = sns.heatmap(file2, ax=axes[1], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap=reversed_cmap, vmin=0, vmax=1,
                cbar_kws={'label': 'FMH dN/dS'}, cbar=False)
axes[1].set_title('10')
axes[1].set_xlabel('')
print(file2.mean())

heatmap3 = sns.heatmap(file3, ax=axes[2], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap=reversed_cmap, vmin=0, vmax=1,
                cbar_kws={'label': 'FMH dN/dS'}, cbar=False)
axes[2].set_title('100')
axes[2].set_xlabel('')
print(file3.mean())


# Create a common color bar
cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])  # Adjust the position and size as needed
cbar = fig.colorbar(heatmap2.get_children()[0], cax=cbar_ax, orientation='vertical', label='FMH dN/dS')


# Adjust layout
#plt.tight_layout()
plt.subplots_adjust(wspace=0.03)
fig.text(0.5, 0.03, 'FMH dN/dS at differing scale factors at ksize=7', ha='center', va='center')

#fig.set_title('ANI and AAI Ground Truth\n10,002 nt at p=0.01')
#fig.figure.savefig(f"{wd}/ANI_and_AAI.png",bbox_inches='tight')
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/fmh_dnds_genome_negative_selection_tests_heatmaps.tiff",bbox_inches='tight') 
