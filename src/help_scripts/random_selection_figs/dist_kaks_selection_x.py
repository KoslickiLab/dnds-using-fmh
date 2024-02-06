import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

p=0.001
selection = 'positive'
wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/kaks_NG'
selection_data = f'{selection}_selection_queries_10002_{p}.axt.kaks'
input_file = pd.read_csv(f'{wd}/{selection_data}',sep='\t')

#input_file = input_file[input_file['Ka/Ks'].notna()] #remove indivisible dNdS ratio
#input_file = input_file[(input_file != 0).all(1)]
#input_file = input_file[(input_file["Ka/Ks"] >= 0) & (input_file["Ka/Ks"] <= 2)]

plt = sns.displot(input_file, x="Ka/Ks",color='red')
plt.set(ylim=(0, 30))
plt.set(xlim=(-0.1, 4))

mean = input_file['Ka/Ks'].mean()
median = input_file['Ka/Ks'].median()
total = len(input_file['Ka/Ks'])

plt.fig.suptitle(f"Distribution of NG dN/dS for {selection} selection (p rate={p}, len=10002)",x=0.5,y=1.03)

ax = plt.facet_axis(0, 0)
ax.text(1.2,22,f'n={total}')
ax.text(1.2,21,f'mean={round(mean,4)}')
ax.text(1.2,20,f'median={round(median,4)}')

ax.axvline(x=1, color='grey', linestyle='--', label='Vertical Line at x=3')
plt.figure.savefig(f'{wd}/{selection}_selection_dist.png',bbox_inches='tight')
