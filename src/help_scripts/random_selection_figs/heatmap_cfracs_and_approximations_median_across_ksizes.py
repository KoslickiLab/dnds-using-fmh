import pandas as pd
import seaborn as sns
from matplotlib import pyplot

#create figure that represents median Cfracs across ksizes for both DNA and Protein
wd='/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/fmh_dnds_sketch_translate'
file = f'{wd}/positive_selection_k5-15/cfracs_and_approximations_median.csv'
clmns=['DNA_Cfrac_median','approx_ANI_median','protein_Cfrac_median','approx_AAI_median','ksize']

new_clmns={"DNA_Cfrac_median":"Positive DNA\nCfrac","protein_Cfrac_median":"Positive AA\nCfrac","approx_ANI_median":"Positive\nANI approx.","approx_AAI_median":"Positive\nAAI approx."}
file1=pd.read_csv(f'{file}',sep=',')[clmns].rename(columns=new_clmns).set_index('ksize')

file = f'{wd}/negative_selection_k5-15/cfracs_and_approximations_median.csv'
new_clmns={"DNA_Cfrac_median":"Negative\nDNA Cfrac","protein_Cfrac_median":"Negative\nAA Cfrac","approx_ANI_median":"Negative\nANI approx.","approx_AAI_median":"Negative\nAAI approx."}
file2=pd.read_csv(f'{file}',sep=',')[clmns].rename(columns=new_clmns).set_index('ksize')

result = pd.concat([file1, file2], axis=1, join="inner")

fig=sns.heatmap(result, annot=True,fmt='.4f',annot_kws={"size": 7})

fig.set_title('Expected Cfrac and Approximated ANI and AAI Medians Across Ksizes\n10,002 nt at p=0.01')
fig.figure.savefig(f"{wd}/cfracs_and_approximations_across_ksizes_medians.png",bbox_inches='tight') 