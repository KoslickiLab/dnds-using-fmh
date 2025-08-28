import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

p=0.1
ksize=7

selection = 'negative' # 'negative' or 'positive'
if selection == 'positive':
    col = 'red'
elif selection == 'negative':
    col='blue'

containment = 'DNA_Cfrac' # 'containment_nt' or 'containment_protein'
if containment == 'containment_nt' or containment == 'DNA_Cfrac':
    sequence = 'DNA'
    k=ksize*3
elif containment == 'containment_protein' or containment == 'AA_Cfrac':
    sequence = 'Protein'
    k=ksize

#wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds/{selection}_selection'
#wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds_sketch_protein/{selection}_selection_redo_sketch_protein_using_faa'
wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds_sketch_protein/{selection}_selection'
selection_data = f'dnds_constant_{ksize}.csv'
input_file = pd.read_csv(f'{wd}/{selection_data}',sep=',')
ref_data = input_file[(input_file["A"] == 'ref_gene')]
other_list = input_file[(input_file["A"] != 'ref_gene')][f'{containment}'].to_list()

mean = ref_data[f'{containment}'].mean()
median = ref_data[f'{containment}'].median()
max = ref_data[f'{containment}'].max()
total = len(ref_data[f'{containment}'])

ref_data_list = ref_data[f'{containment}'].tolist()

# Set the Seaborn style and remove the top and right spines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt = sns.histplot(other_list, label='Pairwise', color='grey', linewidth=0)
plt = sns.histplot(ref_data_list, label='Outgroup', color=f'{col}', linewidth=0)

print(min(ref_data_list))
print(min(other_list))

#print(max(ref_data_list))
#print(max(other_list))

#plt.set_ylim(0, 700)
x_pos=-0.01
#plt.set_xlim(x_pos, 0.1)

font_size=10
add_on=0.05
#plt.text(x_pos+add_on, 400, f'n={total}', fontsize=font_size, color=f'{col}')
#plt.text(x_pos+add_on, 375, f'k={k}', fontsize=font_size, color=f'{col}')
#plt.text(x_pos+add_on, 350, f'mean={round(mean,4)}', fontsize=font_size, color=f'{col}')
#plt.text(x_pos+add_on, 325, f'median={round(median,4)}', fontsize=font_size, color=f'{col}')

plot_title = f"{sequence} Cfrac for {selection} selection (p rate={p}, k={k}, len=10002)"
plt.set_title(plot_title,x=0.5,y=1.03)

legend = plt.legend(loc='upper right')
legend.set_frame_on(False)
for text in legend.get_texts():
    text.set_color('black')



plt.figure.savefig(f'/data/jzr5814/sourmash_dnds_estimation/thesis_figures/containment_distributions/{selection}_{p}_{k}_containment_dist.png',bbox_inches='tight')