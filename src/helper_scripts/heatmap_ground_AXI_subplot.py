import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#create figure that represents median Cfracs across ksizes for both DNA and Protein
wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/ground_truth'
file = f'{wd}/ground_truth_0.01_negative.csv'
clmns=['ANI','AAI']

#new_clmns={"ANI":"Negative\nANI Ground","AAI":"AAI"}
file1=pd.read_csv(f'{file}',sep=',')[clmns]#.rename(columns=new_clmns)

file = f'{wd}/ground_truth_0.01_positive.csv'
#new_clmns={"ANI":"Positive\nANI Ground","AAI":"Positive \nAAI Ground"}
file2=pd.read_csv(f'{file}',sep=',')[clmns]#.rename(columns=new_clmns)

# Create subplots
fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 1]})

# Plot the first heatmap
heatmap1 = sns.heatmap(file1, ax=axes[0], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, vmin=0.94,cmap='Blues',
                cbar_kws={'label': 'Sequence Relatedness'},cbar=False)
axes[0].set_title('Negative')

# Plot the first heatmap
heatmap2 = sns.heatmap(file2, ax=axes[1], annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, vmin=0.94,cmap='Blues',
                cbar_kws={'label': 'Sequence Relatedness'}, cbar=False)
axes[1].set_title('Positive')

# Create a common color bar
cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])  # Adjust the position and size as needed
cbar = fig.colorbar(heatmap2.get_children()[0], cax=cbar_ax, orientation='vertical', label='Sequence Relatedness')


# Adjust layout
#plt.tight_layout()
plt.subplots_adjust(wspace=0.03)
fig.text(0.5, 0.03, 'Ground Truth', ha='center', va='center')

#fig.set_title('ANI and AAI Ground Truth\n10,002 nt at p=0.01')
#fig.figure.savefig(f"{wd}/ANI_and_AAI.png",bbox_inches='tight')
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/ground_truth_ANI_and_AAI_heatmaps.png",bbox_inches='tight') 
