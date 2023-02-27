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

    scatter_png = sns.scatterplot(y='dNdS_ground_truth', x='dNdS_ratio',data=results)
    scatter_png = sns.lmplot(data=results, x="dNdS_ratio", y="dNdS_ground_truth")

    plt.title(f'dNdS_ground_truth (koslick_dNdS) vs. Cfrac dNdS ksize-{ksize}')

    #scatter_png = scatter_png.get_figure()
    scatter_png.figure.savefig(f'{wd}scatter{ksize}.png',bbox_inches='tight')

for k in [5]:
    WD='/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/'
    ground_truth='/data/jzr5814/sourmash_dnds_estimation/tests/data/ground_truth_data/dNdS_ground_truth.csv'
    sourmash_compare_dnds=f'/data/jzr5814/sourmash_dnds_estimation/tests/results/dnds_ground_truth/dNdS_{k}.csv'
    scatter(ground_truth_dnds_results_df=ground_truth,sourmash_compare_dnds_results_df=sourmash_compare_dnds,ksize=k,wd=WD)
