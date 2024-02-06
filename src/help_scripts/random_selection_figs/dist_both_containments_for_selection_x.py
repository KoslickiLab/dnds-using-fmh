import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics

p=0.001
ksize=20

x_pos=1
y_pos=10
substract=1

selection = 'positive'
if selection == 'positive':
    col='red'
elif selection == 'negative':
    col='blue'

wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds'

data_input=f'{wd}/positive_selection/dnds_constant_{ksize}.csv'
data = pd.read_csv(f'{data_input}',sep=',')

if selection == 'positive':
    nt_containment_correct_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] > 1)]['containment_nt'].dropna().tolist()
    prot_containment_correct_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] > 1)]['containment_protein'].dropna().tolist()

    nt_containment_incorrect_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] < 1)]['containment_nt'].dropna().tolist()
    prot_containment_incorrect_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] < 1)]['containment_protein'].dropna().tolist()

elif selection == 'negative':
    nt_containment_correct_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] < 1)]['containment_nt'].dropna().tolist()
    prot_containment_correct_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] < 1)]['containment_protein'].dropna().tolist()

    nt_containment_incorrect_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] > 1)]['containment_nt'].dropna().tolist()
    prot_containment_incorrect_selection = data[(data["A"] == 'ref_gene') & (data['dNdS_ratio_constant'] > 1)]['containment_protein'].dropna().tolist()

plt.hist(nt_containment_correct_selection, alpha=0.5, label='Expected Selection Nucleotide Cfrac', color=f'{col}', linewidth=1.5)
plt.hist(prot_containment_correct_selection, alpha=0.5, label='Expected Negative Selection Protein Cfrac', color=f'{col}', hatch="*", linewidth=1.5)
plt.hist(nt_containment_incorrect_selection, alpha=0.5, label='Unexpected Positive Selection Nucleotide Cfrac', color=f'grey', linewidth=1.5)
plt.hist(prot_containment_incorrect_selection, alpha=0.5, label='Unexpected Positive Selection  Protein Cfrac', color=f'grey', hatch="*", linewidth=1.5)

ax = plt.gca()                        
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

if nt_containment_correct_selection:
    nt_containment_correct_selection_median = statistics.median(sorted(nt_containment_correct_selection))
    plt.text(x_pos,y_pos,f'nt Cfrac expected selection median={round(nt_containment_correct_selection_median, 4)}')
else:
    plt.text(x_pos,y_pos,f'nt Cfrac expected selection median=nan')

if nt_containment_incorrect_selection:
    nt_containment_incorrect_selection_median = statistics.median(sorted(nt_containment_incorrect_selection))
    plt.text(x_pos,y_pos-substract,f'nt Cfrac unexpected selection median={round(nt_containment_incorrect_selection_median, 4)}')
else:
    plt.text(x_pos,y_pos-substract,f'nt Cfrac unexpected selection median=nan')

if prot_containment_correct_selection:
    prot_containment_correct_selection_median = statistics.median(sorted(prot_containment_correct_selection))
    plt.text(x_pos,y_pos-substract-1,f'protein Cfrac expected selection median={round(prot_containment_correct_selection_median, 4)}')
else:
    plt.text(x_pos,y_pos-substract-1,f'protein Cfrac expected selection median=nan')

if prot_containment_incorrect_selection:
    prot_containment_incorrect_selection_median = statistics.median(sorted(prot_containment_incorrect_selection))
    plt.text(x_pos,y_pos-substract-2,f'protein Cfrac unexpected selection median={round(prot_containment_incorrect_selection_median, 4)}')
else:
    plt.text(x_pos,y_pos-substract-2,f'protein Cfrac unexpected selection median=nan')


plot_title = f"Distribution of nt and protein Cfrac for {selection} selection \np rate={p}, len=10002, k={ksize}"
plt.suptitle(plot_title,x=0.5,y=0.95)

plt.legend(bbox_to_anchor=(1, 1))


plt.savefig(f'{wd}/{ksize}_correct_and_incorrect_{selection}_for_both_containments_dist.png',bbox_inches='tight')
