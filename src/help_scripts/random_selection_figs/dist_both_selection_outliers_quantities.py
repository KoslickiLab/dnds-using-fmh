import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics

p=0.001
selection = 'negative'
if selection == 'positive':
    color = 'red'
elif selection == 'negative':
    color='blue'

ksize=7
wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds_not_using_sourmash_translate'

data_positive=f'{wd}/positive_selection/dnds_constant_{ksize}.csv'
positive = pd.read_csv(f'{data_positive}',sep=',')
positive = positive[(positive["A"] == 'ref_gene') & (positive["dNdS_ratio_constant"] < 0) | (positive["dNdS_ratio_constant"] > 2)][f'dNdS_ratio_constant'].dropna().tolist()

data_negative=f'{wd}/negative_selection/dnds_constant_{ksize}.csv'
negative = pd.read_csv(f'{data_negative}',sep=',')
negative = negative[(negative["A"] == 'ref_gene') & (negative["dNdS_ratio_constant"] < 0) | (negative["dNdS_ratio_constant"] > 2)][f'dNdS_ratio_constant'].dropna().tolist()

# Set the Seaborn style and remove the top and right spines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt = sns.histplot(positive, label='Positive', color='red', linewidth=0)
plt = sns.histplot(negative, label='Negative', color='blue', linewidth=0)

#plt.set(ylim=(0, 30))
#plt.set(xlim=(-0.1, 4))

#mean = positive['dNdS_ratio_constant'].mean()

#total = len(positive['dNdS_ratio_constant'])

plot_title = f"Distribution of FMH dNdS Outliers between positive and negative selection (p rate={p}, len=10002, k={ksize})"
plt.set_title(plot_title,x=0.5,y=1.03)

#ax = plt.facet_axis(0, 0)
#ax.text(1.2,29,f'k={ksize}')
#ax.text(1.2,28,f'n={total}')
#ax.text(1.2,27,f'mean={round(mean,4)}')
#ax.text(1.2,26,f'median={round(median,4)}')

plt.text(1.1,25,f'k={ksize}')

positive_median = statistics.median(sorted(positive))
plt.text(1.1,24,f'positive selection median={round(positive_median, 4)}')

negative_median = statistics.median(sorted(negative))
plt.text(1.1,23,f'negative selection median={round(negative_median, 4)}')

plt.axvline(x=1, color='grey', linestyle='--', label='Vertical Line at x=3')
plt.figure.savefig(f'{wd}/{ksize}_dnds_selection_outlier_dist.png',bbox_inches='tight')
