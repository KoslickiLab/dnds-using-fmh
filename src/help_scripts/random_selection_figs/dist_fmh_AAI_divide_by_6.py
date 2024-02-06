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

ksize=7
wd = f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.001/divide_by_6'
data=f'{wd}/AAI_divide_by_6.csv'
df = pd.read_csv(f'{data}',sep=',')

df["AAI/6"] = df["AAI"]/6
plt = sns.displot(df, x="AAI/6",color=color, linewidth=0)
#plt.set(ylim=(0, 30))
#plt.set(xlim=(-0.1, 4))


mean = df["AAI/6"].mean()
median = df["AAI/6"].median()
total = len(df["AAI/6"])

plt.fig.suptitle(f"AAI/6 sourmash translate\np rate={p}, len=10002, k_aa={ksize}",x=0.5,y=1.03)

ax = plt.facet_axis(0, 0)
#ax.text(-0.352,29,f'k={ksize}')
#ax.text(-0.352,28,f'n={total}')
#x.text(-0.352,27,f'median={round(median,4)}')
#ax.text(1.2,26,f'median={round(median,4)}')


#ax.axvline(x=1, color='grey', linestyle='--', label='Vertical Line at x=3')
plt.set_xticklabels(rotation=40)

plt.figure.savefig(f'{wd}/{selection}_{ksize}_selection_AAI_divide_by_6_dist3.png',bbox_inches='tight')
