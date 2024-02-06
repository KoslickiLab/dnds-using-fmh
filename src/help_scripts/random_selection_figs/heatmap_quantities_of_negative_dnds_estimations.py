import pandas as pd
import seaborn as sns
from matplotlib import pyplot

#create figure that represents median Cfracs across ksizes for both DNA and Protein
wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/fmh_dnds_sketch_protein/positive_selection'
file = f'{wd}/dnds_constant_all_fixed_columns.csv'
clmns=['A','B','ksize','DNA_Cfrac','AA_Cfrac','PdN','PdS','dN/dS']
k=7
quantity='dN/dS'

file1=pd.read_csv(f'{file}',sep=',')[clmns]
print(file1)
file1["dN/dS"] = pd.to_numeric(file1["dN/dS"])
file1["ksize"] = pd.to_numeric(file1["ksize"])

file1['approx_ANI'] = file1['DNA_Cfrac']**(1/(3*file1['ksize']))
file1['approx_AAI'] = file1['AA_Cfrac']**(1/file1['ksize'])

file1=file1[(file1['A']=='ref_gene')&(file1['ksize']==k)&(file1['dN/dS']<0)][['B','approx_ANI','approx_AAI']].set_index('B')
#file1=file1[(file1['dN/dS']<0)]
print(file1)




fig=sns.heatmap(file1, annot=False,fmt='.4f',annot_kws={"size": 7},yticklabels=False)

fig.set_title(f'ANI and AAI Approximations\nPositive, 10,002 nt, p=0.001, k={k}')
fig.figure.savefig(f"{wd}/ANI_and_AAI_approx_{k}.png",bbox_inches='tight') 