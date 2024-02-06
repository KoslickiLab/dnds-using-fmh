import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

p=0.001
selection = 'negative'

if selection == 'positive':
    color = 'red'
elif selection == 'negative':
    color='blue'

ksize=5
wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds_not_using_sourmash_translate/{selection}_selection/'
data_positive=f'{wd}dnds_constant_{ksize}.csv'
positive = pd.read_csv(f'{data_positive}',sep=',')

#positive = positive.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
#positive = positive[(positive != 0).all(1)]
positive = positive[(positive["dNdS_ratio_constant"] >= 0) & (positive["dNdS_ratio_constant"] <= 2)]
positive = positive[(positive["A"] == 'ref_gene')]

plt = sns.displot(positive, x="dNdS_ratio_constant",color=color)
#plt.set(ylim=(0, 30))
#plt.set(xlim=(-0.1, 4))


mean = positive['dNdS_ratio_constant'].mean()
median = positive['dNdS_ratio_constant'].median()
total = len(positive['dNdS_ratio_constant'])

plt.fig.suptitle(f"Distribution of FMH dN/dS for {selection} selection p rate={p}, len=10002, k={ksize}",x=0.5,y=1.03)

ax = plt.facet_axis(0, 0)
ax.text(1.2,29,f'k={ksize}')
ax.text(1.2,28,f'n={total}')
ax.text(1.2,27,f'mean={round(mean,4)}')
#ax.text(1.2,26,f'median={round(median,4)}')


ax.axvline(x=1, color='grey', linestyle='--', label='Vertical Line at x=3')
plt.figure.savefig(f'{wd}{selection}_{ksize}_selection_dist.png',bbox_inches='tight')
