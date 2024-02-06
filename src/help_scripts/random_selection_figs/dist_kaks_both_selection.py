import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics

p=0.01

wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/kaks_NG'

data_positive = f'positive_selection_queries_10002_{p}.axt.kaks'
positive = pd.read_csv(f'{wd}/{data_positive}',sep='\t')
positive = positive['Ka/Ks'].dropna().tolist()
print(positive)

data_negative= f'negative_selection_queries_10002_{p}.axt.kaks'
negative = pd.read_csv(f'{wd}/{data_negative}',sep='\t')
negative = negative['Ka/Ks'].dropna().tolist()
print(negative)

# Set the Seaborn style and remove the top and right spines
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt = sns.histplot(positive, label='Positive', color='red', linewidth=0)
plt = sns.histplot(negative, label='Negative', color='blue', linewidth=0)

#plt.set(ylim=(0, 30))

plot_title = f"NG86-Based dN/dS of Selectively Mutated Sequences\np-rate={p}, len=10002"
plt.set_title(plot_title,x=0.5,y=1.03)

#plt.text(1.2,29,f'k={ksize}')

positive_median = statistics.median(sorted(positive))

#plt.text(6,20,f'positive selection median={round(positive_median, 4)}')

negative_median = statistics.median(sorted(negative))
#plt.text(6,19,f'negative selection median={round(negative_median, 4)}')
ax.set_xlabel('NG86 dN/dS')
ax.legend()
plt.axvline(x=1, color='grey', linestyle='--', label='Vertical Line at x=3')
plt.figure.savefig(f'{wd}/kaks_selection_dist_thesis.png',bbox_inches='tight')
