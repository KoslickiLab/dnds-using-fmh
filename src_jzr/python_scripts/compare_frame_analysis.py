import pandas as pd
import csv
import pickle, os, glob, subprocess

def frame(x):
    return(x[-1:])

wd = '/data/jzr5814/data/frame_analysis_using_bash_script/compare_by_multiple_jobs/test1/'

ksizes = [7,14,21,28,35,42,49,56,63,70]

for ksize in ksizes:
    df = pd.read_csv(wd+'compare'+str(ksize)+'.csv', sep=',')
    transpose_df=df.T.iloc[:,0].reset_index().rename(columns={'index':'query',0:'containment'}) #extract column where all queries are compared to ref
    transpose_df['ksize'] = ksize
#    transpose_df['frame'] = transpose_df["query"][-1:]

    transpose_df["frame"] = transpose_df["query"].astype(str).apply(frame)


    transpose_df.to_csv(wd+'compare_containment'+str(ksize)+'.csv',sep=',') 

cmd = f"awk '(NR == 1) || (FNR > 1)' {wd}compare_containment*.csv | grep -v '../../data/uniprotkb.fasta' > {wd}containment.csv"
subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)


