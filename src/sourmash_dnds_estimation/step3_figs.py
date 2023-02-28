import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pylab import rcParams

def scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,ksize,wd):
    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    results = pd.concat([sourmash_compare_dnds_results_df['dNdS_ratio'],expected_dnds_results_df['dNdS_ground_truth']], axis=1).reindex(expected_dnds_results_df.index)
    #minimum = results.min().min()
    #maximum = results.max().max()

    scatter_png = sns.scatterplot(x='dNdS_ratio',y='dNdS_ground_truth',data=results)
    scatter_png = sns.lmplot(data=results, x="dNdS_ratio", y="dNdS_ground_truth")

    plt.title(f'dNdS_ground_truth (koslick_dNdS) vs. Cfrac dNdS ksize={ksize}')
    plt.xlabel('containment_dNdS')
    #scatter_png = scatter_png.get_figure()
    scatter_png.figure.savefig(f'{wd}scatter{ksize}.png',bbox_inches='tight')

for k in [2,3,4,5,6,7,8,9,10,11,12,13,14,15,20]:
    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/'
    ground_truth='/data/jzr5814/sourmash_dnds_estimation/tests/data/ground_truth_data/dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/dNdS_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD)

#for k in [2,3,4,5,6,7,8,9,10,15,20]:
#    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/'
#    ground_truth='/data/jzr5814/sourmash_dnds_estimation/tests/data/ground_truth_data/dNdS_ground_truth.csv'
#    sourmash_compare_dnds=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/dNdS_{k}.csv'
#    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD)

def multiplot_scatter(ground_truth_dnds_results_df,sourmash_compare_dnds_results_df,ksize):
    expected_dnds_results_df = pd.read_csv(ground_truth_dnds_results_df)
    sourmash_compare_dnds_results_df = pd.read_csv(sourmash_compare_dnds_results_df)

    results = pd.concat([sourmash_compare_dnds_results_df['dNdS_ratio'],expected_dnds_results_df['dNdS_ground_truth']], axis=1).reindex(expected_dnds_results_df.index)
    #minimum = results.min().min()
    #maximum = results.max().max()

    scatter_png = sns.scatterplot(x='dNdS_ratio',y='dNdS_ground_truth',data=results)
    scatter_png = sns.lmplot(data=results, x="dNdS_ratio", y="dNdS_ground_truth")

    plt.title(f'dNdS_ground_truth (koslick_dNdS) vs. Cfrac dNdS ksize={ksize}')
    plt.xlabel('containment_dNdS')

    return(scatter_png)

fig, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2)
fig.suptitle('A tale of 2 subplots')

ax1.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=5)
ax2.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=6)
ax3.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=7)
ax4.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=8)
ax5.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=9)
ax6.multiplot_scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=10)


WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/'
plt.figure.savefig(f'{WD}test.png')



