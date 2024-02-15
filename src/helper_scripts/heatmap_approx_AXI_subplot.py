import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#create figure that represents median Cfracs across ksizes for both DNA and Protein
clmns=['A','B','ksize','DNA_Cfrac','AA_Cfrac','PdN','PdS','dN/dS']
quantity='dN/dS'
k=7
p=0.01


wd1=f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds_sketch_protein/positive_selection_redo_sketch_protein_using_faa'
file1 = f'{wd1}/all_dnds_constant.csv'
file1=pd.read_csv(f'{file1}',sep=',')[clmns]
file1["ksize"] = pd.to_numeric(file1["ksize"])
file1['approx_ANI'] = file1['DNA_Cfrac']**(1/(3*file1['ksize']))
file1['approx_AAI'] = file1['AA_Cfrac']**(1/file1['ksize'])
file1=file1[(file1['A']=='ref_gene')&(file1['ksize']==k)][['approx_ANI','approx_AAI']].rename(columns={'approx_ANI':'ANI','approx_AAI':'AAI'})

wd2=f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds_sketch_protein/negative_selection_redo_sketch_protein_using_faa'
file2 = f'{wd1}/all_dnds_constant.csv'
file2=pd.read_csv(f'{file2}',sep=',')[clmns]
print(file2)
file2["ksize"] = pd.to_numeric(file2["ksize"])
file2['approx_ANI'] = file2['DNA_Cfrac']**(1/(3*file2['ksize']))
file2['approx_AAI'] = file2['AA_Cfrac']**(1/file2['ksize'])
file2=file2[(file2['A']=='ref_gene')&(file2['ksize']==k)][['approx_ANI','approx_AAI']].rename(columns={'approx_ANI':'ANI','approx_AAI':'AAI'})

# Create subplots
fig, axes = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 1]})

# Plot the first heatmap
heatmap1 = sns.heatmap(file2, ax=axes[0], annot=False,fmt='.4f',annot_kws={"size": 7},yticklabels=False,cmap='Blues',vmin=0.94,cbar_kws={'label': 'Sequence Relatedness'},cbar=False)
axes[0].set_title('Negative')

# Plot the first heatmap
heatmap2 = sns.heatmap(file1, ax=axes[1], annot=False,fmt='.4f',annot_kws={"size": 7},yticklabels=False, cmap='Blues',vmin=0.94,cbar_kws={'label': 'Sequence Relatedness'}, cbar=False)
axes[1].set_title('Positive')

# Create a common color bar
cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])  # Adjust the position and size as needed
cbar = fig.colorbar(heatmap2.get_children()[0], cax=cbar_ax, orientation='vertical', label='Sequence Relatedness')

fig.text(0.5, 0.03, f'Approximations from FMH Containment {p},k={k}', ha='center', va='center')
# Adjust layout
#plt.tight_layout()
plt.subplots_adjust(wspace=0.03)

#fig.set_title('ANI and AAI Ground Truth\n10,002 nt at p=0.01')
#fig.figure.savefig(f"{wd}/ANI_and_AAI.png",bbox_inches='tight')
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/approx_ANI_and_AAI_{p}_heatmaps.png",bbox_inches='tight') 
