import pandas as pd
import seaborn as sns
from matplotlib import pyplot

#create figure that represents median Cfracs across ksizes for both DNA and Protein
wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/ground_truth'
file = f'{wd}/ground_truth_0.01_negative.csv'
clmns=['ANI','AAI']

new_clmns={"ANI":"Negative\nANI Ground","AAI":"Negative\nAAI Ground"}
file1=pd.read_csv(f'{file}',sep=',')[clmns].rename(columns=new_clmns)

file = f'{wd}/ground_truth_0.01_positive.csv'
new_clmns={"ANI":"Positive\nANI Ground","AAI":"Positive \nAAI Ground"}
file2=pd.read_csv(f'{file}',sep=',')[clmns].rename(columns=new_clmns)

result = pd.concat([file1, file2], axis=1, join="inner")
print(result)
#sns.color_palette("Blues", as_cmap=True)
fig=sns.heatmap(result, annot=False,fmt='.4f',annot_kws={"size": 7}
                ,yticklabels=False, vmin=0.94,cmap='Blues',
                cbar_kws={'label': 'Sequence Relatedness'})

#fig.set_title('ANI and AAI Ground Truth\n10,002 nt at p=0.01')
#fig.figure.savefig(f"{wd}/ANI_and_AAI.png",bbox_inches='tight')
fig.figure.savefig(f"/data/jzr5814/sourmash_dnds_estimation/thesis_figures/ground_truth_ANI_and_AAI_heatmaps.png",bbox_inches='tight') 