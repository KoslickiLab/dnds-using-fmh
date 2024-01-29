import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr

p = 0.01
select ='positive'
len=10002

AAI_file = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/ground_truth/ground_truth_0.01_positive.csv",sep=",",header=0)
dS_file = pd.read_csv(f"/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/kaks_NG/positive_selection_queries_10002_0.01.axt.kaks",sep="\t",header=0)


print(AAI_file)

print(dS_file)
#create figure    
plt.scatter(AAI_file["AAI"],dS_file["Ks"],color='red')
plt.title(f'AAI vs Synonymous Mutation Rate\n{select} {p} len={len}')
plt.xlabel(f'AAI')
plt.ylabel(f'dS')
#plt.plot(min(ANI_file["ANI"]),max(dS_file["Ks"]),linestyle="--")
plt.savefig(f'/data/jzr5814/sourmash_dnds_estimation/tests/test/create_sequence_using_NG_assumption/0.01/scatter_AAI_to_Ks_{select}_{p}.png')
