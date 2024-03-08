import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#create figure that represents median Cfracs across ksizes for both DNA and Protein
clmns=['dN/dS']

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k7_1scale_cores1000/fmh_omega_7.csv'
file1=pd.read_csv(f'{file}',sep=',')[clmns]
file1=file1[file1 > 0]

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k7_10scale_cores1000/fmh_omega_7.csv'
file2=pd.read_csv(f'{file}',sep=',')[clmns]
file2=file2[file2 > 0]

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k7_100scale_cores1000/fmh_omega_7.csv'
file3=pd.read_csv(f'{file}',sep=',')[clmns]
file3=file3[file3 > 0]

file='/data/jzr5814/sourmash_dnds_estimation/tests/test/genome_selection/positive_k7_1000scale_cores1000/fmh_omega_7.csv'
file4=pd.read_csv(f'{file}',sep=',')[clmns]
file4=file4[file4 > 0]

max_value = max(max(file1['dN/dS']),max(file2['dN/dS']),max(file3['dN/dS']),max(file4['dN/dS']))
#if max_value > 10:
#    file1=file1[(file1 <= 10) & (file1 > 0)]
#    file2=file2[(file2 <= 10) & (file2 > 0)]
#    file3=file3[(file3 <= 10) & (file3 > 0)]
#    file4=file4[(file4 <= 10) & (file4 > 0)]
#    max_value=10

# Create subplots
fig, axes = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1, 1, 1, 1]})

# Get the colormap and reverse it
cmap = plt.get_cmap('vlag')

# Plot the first heatmap
heatmap1 = sns.heatmap(file1, ax=axes[0], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap=cmap, vmin=0, vmax=max_value, center=1,
                cbar_kws={'label': 'FMH dN/dS'},cbar=False)
axes[0].set_title('1')
axes[0].set_xticks([])
print(file1.describe())

heatmap2 = sns.heatmap(file2, ax=axes[1], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap=cmap, vmin=0, vmax=max_value, center=1,
                cbar_kws={'label': 'FMH dN/dS'}, cbar=False)
axes[1].set_title('10')
axes[1].set_xticks([])
print(file2.describe())

heatmap3 = sns.heatmap(file3, ax=axes[2], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap=cmap, vmin=0, vmax=max_value, center=1,
                cbar_kws={'label': 'FMH dN/dS'}, cbar=False)
axes[2].set_title('100')
axes[2].set_xticks([])
print(file3.describe())

heatmap4 = sns.heatmap(file4, ax=axes[3], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, cmap=cmap, vmin=0, vmax=max_value, center=1,
                cbar_kws={'label': 'FMH dN/dS'}, cbar=False)
axes[3].set_title('1000')
axes[3].set_xticks([])
print(file4.describe())


# Create a common color bar
cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])  # Adjust the position and size as needed
cbar = fig.colorbar(heatmap2.get_children()[0], cax=cbar_ax, orientation='vertical', label='FMH dN/dS')


# Adjust layout
plt.subplots_adjust(wspace=0.03)
plt.suptitle('positive selection simulation on E. coli protein rep genome, p=0.01,ksize=7')
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/fmh_dnds_genome_positive_selection_tests_heatmaps.tiff",bbox_inches='tight') 
