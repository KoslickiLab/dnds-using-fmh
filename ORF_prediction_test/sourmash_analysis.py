import pickle
import os
import glob
import pandas as pd

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

smallk = pd.read_csv("../data/results_kmers.csv",sep=",",header=None,engine="python",on_bad_lines="skip",quoting=3)

def frame(x):
    return(x[-7:])
def sample(x):
    return(x[:15])

match_name=1
max_containment=0
ksize=2

smallk["frame"] = smallk[match_name].apply(frame)
smallk = smallk[smallk['frame'].str.startswith('frame_')]
smallk['sample'] = smallk[match_name].apply(sample)
smallk[ksize] = pd.to_numeric(smallk[ksize])
smallk = smallk[["frame","sample",max_containment,ksize]].rename(columns={0: 'max_containment', 2: 'ksize'})
#smallk = smallk.sort_values(by=['sample', max_containment],ascending=False).rename(columns={0: 'max_containment', 2: 'ksize'})
smallk = smallk.sort_values(by=['sample', "ksize", "max_containment"],ascending=False).rename(columns={0: 'max_containment', 2: 'ksize'})
#print(smallk.groupby(["sample","ksize"]).head(1).reset_index(drop=True))
print(smallk)
smallk.to_csv('../data/containment.csv')
#print(smallk.groupby("sample").head(2).reset_index(drop=True))
#print(smallk.groupby(["ksize","max_containment"]).head(2))
#print(smallk.shape)

temp_sample = ''
ignore_sample=''
dictionary={}
for index, row in smallk.iterrows():
    if temp_sample!=row['sample']:
        temp_sample = row['sample']
        if row["ksize"] not in dictionary:
            if row['frame'] == 'frame_1':
                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                dictionary[row['ksize']]['highest'].append(row['max_containment'])
            else:
                if ignore_sample != temp_sample:
                    dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
                    ignore_sample=temp_sample
        elif row["ksize"] in dictionary:
            if row['frame'] == 'frame_1':
                dictionary[row['ksize']]['highest'].append(row['max_containment'])
            else:
                if ignore_sample != temp_sample:
                    dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
                    ignore_sample=temp_sample
    elif temp_sample==row['sample']:
        if row["ksize"] not in dictionary:
            if row['frame'] == 'frame_1':
                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                dictionary[row['ksize']]['highest'].append(row['max_containment'])
            else:
                dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
        elif row["ksize"] in dictionary:
            if row['frame'] == 'frame_1':
                dictionary[row['ksize']]['highest'].append(row['max_containment'])
            else:
                if ignore_sample != temp_sample:
                    dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
                    ignore_sample=temp_sample
print(dictionary)

with open('dictionary.pickle', 'wb') as handle:
    pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)



