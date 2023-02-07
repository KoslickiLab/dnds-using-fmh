import pandas as pd
import csv

wd = '/data/jzr5814/data/frame_analysis_using_bash_script/compare_by_multiple_jobs/test1/'
ksizes = [7,14,21,28,35,42,49,56,63,70]

for ksize in ksizes:
    df = pd.read_csv(wd+'compare'+str(ksize)+'.csv', sep=',')
    transpose_df=df.T.iloc[:,0] #extract column where all queries are compared to ref
    transpose_df.to_csv(wd+'compare_containment'+str(ksize)+'.csv',header=False,sep=',') 
