import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

p=0.001
selection = 'positive'

if selection == 'positive':
    color = 'red'
elif selection == 'negative':
    color='blue'

ksize=5
wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/{p}/fmh_dnds_not_using_sourmash_translate/{selection}_selection/'
data_positive=f'{wd}dnds_constant_{ksize}.csv'
positive = pd.read_csv(f'{data_positive}',sep=',')
positive = positive[(positive["A"] == 'ref_gene')]
#positive = positive.replace([np.inf, -np.inf], np.nan).dropna() #remove indivisible dNdS ratio
#positive = positive[(positive != 0).all(1)]

nan_count = positive['dNdS_ratio_constant'].isna().sum()

positive = positive[(positive["dNdS_ratio_constant"] >= 0) & (positive["dNdS_ratio_constant"] <= 2)]
#positive = positive[(positive["dNdS_ratio_constant"] < 0)]
#positive = positive[(positive["dNdS_ratio_constant"] > 2)]

plt = sns.displot(positive, x="dNdS_ratio_constant",color=color, linewidth=0)
#plt.set(ylim=(0, 10))
#plt.set(xlim=(-50, 100))


mean = positive['dNdS_ratio_constant'].mean()
median = positive['dNdS_ratio_constant'].median()
total = len(positive['dNdS_ratio_constant'])
print(total)

plt.fig.suptitle(f"FMH dN/dS estimations outliers included\n{selection} selection, p rate={p}, len=10002, k={ksize}\n",x=0.5,y=1.03)
#plt.fig.suptitle(f"Distribution of FMH dN/dS below zero for {selection} selection\np rate={p}, len=10002, k={ksize}",x=0.5,y=1.03)

ax = plt.facet_axis(0, 0)
#ax.text(1000,5,f'k={ksize}',color='white')
ax.text(0.75,30,f'n={total}',color='black')
ax.text(0.75,28,f'mean={round(mean,4)}',color='black')
ax.text(0.75,26,f'median={round(median,4)}',color='black')


#ax.axvline(x=0, color='grey', linestyle='--', label='Vertical Line at x=3')
#plt.figure.savefig(f'{wd}{selection}_{ksize}_dnds_selection_outliers_excluded_dist.png',bbox_inches='tight')
