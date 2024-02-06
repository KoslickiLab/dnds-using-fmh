import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics

p=0.001
ksize=15

containment = 'containment_protein' # 'containment_nt' or 'containment_protein'
if containment == 'containment_nt':
    sequence = 'nucleotide'
    k=ksize*3
elif containment == 'containment_protein':
    sequence = 'protein'
    k=ksize

wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds'

data_positive=f'{wd}/positive_selection/dnds_constant_{ksize}.csv'
positive = pd.read_csv(f'{data_positive}',sep=',')
positive_correct = positive[(positive["A"] == 'ref_gene') & (positive['dNdS_ratio_constant'] > 1)][f'{containment}'].dropna().tolist()
positive_incorrect = positive[(positive["A"] == 'ref_gene') & (positive['dNdS_ratio_constant'] < 1)][f'{containment}'].dropna().tolist()

data_negative=f'{wd}/negative_selection/dnds_constant_{ksize}.csv'
negative = pd.read_csv(f'{data_negative}',sep=',')
negative_correct = negative[(negative["A"] == 'ref_gene') & (negative['dNdS_ratio_constant'] < 1)][f'{containment}'].dropna().tolist()
negative_incorrect = negative[(negative["A"] == 'ref_gene') & (negative['dNdS_ratio_constant'] > 1)][f'{containment}'].dropna().tolist()

plt.hist(positive_correct, alpha=0.5, label='Expected Positive Selection', color='red', linewidth=1.5)
plt.hist(negative_correct, alpha=0.5, label='Expected Negative Selection', color='blue', linewidth=1.5)
plt.hist(positive_incorrect, alpha=0.5, label='Unexpected Positive Selection', color='red',hatch="/", linewidth=1.5)
plt.hist(negative_incorrect, alpha=0.5, label='Unexpected Positive Selection', color='blue', hatch="/", linewidth=1.5)

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

x_pos=1
if positive_correct:
    positive_correct_median = statistics.median(sorted(positive_correct))
    plt.text(x_pos,9,f'expected positive selection median={round(positive_correct_median, 4)}')
else:
    plt.text(x_pos,9,f'expected positive selection median=nan')

if positive_incorrect:
    positive_incorrect_median = statistics.median(sorted(positive_incorrect))
    plt.text(x_pos,8,f'unexpected positive selection median={round(positive_incorrect_median, 4)}')
else:
    plt.text(x_pos,8,f'unexpected positive selection median=nan')

if negative_correct:
    negative_correct_median = statistics.median(sorted(negative_correct))
    plt.text(x_pos,7,f'expected negative selection median={round(negative_correct_median, 4)}')
else:
    plt.text(x_pos,7,f'expected negative selection median=nan')

if negative_incorrect:
    negative_incorrect_median = statistics.median(sorted(negative_incorrect))
    plt.text(x_pos,6,f'unexpected negative selection median={round(negative_incorrect_median, 4)}')
else:
    plt.text(x_pos,6,f'unexpected negative selection median=nan')

plot_title = f"Distribution of {sequence} containments for positive and negative selection \np rate={p}, len=10002, k={k}"
plt.suptitle(plot_title,x=0.5,y=0.95)

plt.legend(bbox_to_anchor=(1, 1))

"""
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
plt.text(x_pos+add_on, 400, f'n={total}', fontsize=font_size, color=f'{col}')
plt.text(x_pos+add_on, 375, f'k={k}', fontsize=font_size, color=f'{col}')
plt.text(x_pos+add_on, 350, f'mean={round(mean,4)}', fontsize=font_size, color=f'{col}')
plt.text(x_pos+add_on, 325, f'median={round(median,4)}', fontsize=font_size, color=f'{col}')

plot_title = f"Distribution of {sequence} containment for {selection} selection (p rate={p}, len=10002)"
plt.set_title(plot_title,x=0.5,y=1.03)

legend = plt.legend(loc='upper right')
legend.set_frame_on(False)
for text in legend.get_texts():
    text.set_color('black')


"""
plt.savefig(f'{wd}/{containment}_{k}_correct_and_incorrect_selection_dist.png',bbox_inches='tight')
