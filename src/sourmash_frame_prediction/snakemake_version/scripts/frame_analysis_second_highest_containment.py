import pandas as pd
import csv
import pickle, os, glob, subprocess

wd='/data/jzr5814/data/frame_analysis_using_bash_script/compare_by_multiple_jobs/test1/'

for k in [7,14,21,28,35,42,49,56,63,70]:
    df=pd.read_csv(wd+'compare_containment'+str(k)+'.csv',sep=',')[['query','ksize','containment','frame']].iloc[1:]
    df=df[df['frame']!='frame_1']
    df['gene']=df['query'].apply(lambda x: x[:-2])

    dict={}
    for index,row in df.iterrows():
        if row['gene'] not in dict:
                dict[row['gene']]={}
                dict[row['gene']]['ksize']=[]
                dict[row['gene']]['ksize'].append(row['ksize'])
                dict[row['gene']]['containment']=[]
                dict[row['gene']]['containment'].append(row['containment'])
        elif row['gene'] in dict:
                dict[row['gene']]['ksize'].append(row['ksize'])
                dict[row['gene']]['containment'].append(row['containment'])

    new_dict={}
    new_dict['gene']=[]
    new_dict['ksize']=[]
    new_dict['containment']=[]
    for key in dict:
        CI = max(dict[key]['containment'])
        ksize = dict[key]['ksize'][dict[key]['containment'].index(max(dict[key]['containment']))]
        if key not in new_dict:
            new_dict['gene'].append(key)
            new_dict['ksize'].append(ksize)
            new_dict['containment'].append(CI)

    new_df=pd.DataFrame.from_dict(new_dict)
    new_df.describe()

    pd.DataFrame.from_dict(new_dict).rename_axis('index').to_csv(wd+'second_highest_containment'+str(k)+'.csv')

cmd1 = f"awk '(NR == 1) || (FNR > 1)' {wd}second_highest_containment*.csv > {wd}second_highest_containment.csv"
subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)
