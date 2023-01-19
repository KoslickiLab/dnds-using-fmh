import pickle, os, glob, subprocess
import pandas as pd

def produce_containment_csv(wd):
    cmd1 = f"awk '(NR == 1) || (FNR > 1)' {wd}/prefetch_res_*.csv | cut -d ',' -f 3,7,14,15 | tr '_' ',' | sed 's/max,containment/max_containment/' > {wd}results_kmers.csv"
    subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

def frame(x):
#    print(x)
    return('frame_'+x)

def prep(data,output="containment.csv"):
    kmers_data = pd.read_csv(data,sep=",",converters={"name":str}).rename(columns = {'match':'query','name':'frame'}).sort_values(by=["query", "ksize", "max_containment"],ascending=False)
    #match_name=1
    #max_containment=0
    kmers_data["frame"] = kmers_data["frame"].astype(str).apply(frame)

    kmers_data.to_csv(output)

def analysis(data,output="CIdict.pickle"):
    temp_sample = ''
    ignore_sample=''
    dictionary={}
    kmers_data = pd.read_csv(data,sep=",")
    for index, row in kmers_data.iterrows():
        if temp_sample!=row["query"]:
            temp_sample = row["query"]
            if row["ksize"] not in dictionary:
                if row['frame'] == 'frame_01':
                    dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
                else:
                    if ignore_sample != temp_sample:
                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
                        ignore_sample=temp_sample
            elif row["ksize"] in dictionary:
                if row['frame'] == 'frame_01':
                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
                else:
                    if ignore_sample != temp_sample:
                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
                        ignore_sample=temp_sample
        elif temp_sample==row["query"]:
            if row["ksize"] not in dictionary:
                if row['frame'] == 'frame_01':
                    dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                else:
                    dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
            elif row["ksize"] in dictionary:
                if row['frame'] == 'frame_01':
                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
                else:
                    if ignore_sample != temp_sample:
                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
                        ignore_sample=temp_sample

    with open(output, 'wb') as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)



