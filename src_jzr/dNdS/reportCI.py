import pickle, os, glob, subprocess
import pandas as pd

def produce_containment_csv(wd):
    cmd1 = f"awk '(NR == 1) || (FNR > 1)' {wd}prefetch_res_*.csv | cut -d ',' -f 3,7,14,15 | tr '_' ',' | sed 's/max,containment/max_containment/' > {wd}results_kmers.csv"
    subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

def frame(x):
#    print(x)
    return('frame_'+x)

def prep(data,output="containment.csv"):
    kmers_data = pd.read_csv(data,sep=",",converters={"name":str}).rename(columns = {'match':'query','name':'frame'}).sort_values(by=["query", "ksize", "max_containment"],ascending=False)
    #match_name=1
    #max_containment=0
    kmers_data["frame"] = kmers_data["frame"].astype(str).apply(frame)

    kmers_data.to_csv(output,encoding='utf-8')

def extract_frame1(data="containment.csv",output="frame1_containment.csv"):
    #Produce csv file of frame 1 containmemt indexes
    frame1_data = pd.read_csv(data,sep=",").reset_index()
    frame1_data = frame1_data[frame1_data['frame'] == 'frame_1']
    frame1_data[['index','max_containment','frame','ksize']].to_csv(output,encoding='utf-8',index=False)

def extract_frameX(data="containment.csv",output="frameX_containment.csv"):
    #Produce csv file of excluded frame 1 containmemt indexes
    framex_data = pd.read_csv(data,sep=",").reset_index()
    framex_data = framex_data[framex_data['frame'] != 'frame_1']
    framex_data[['index','max_containment','frame','ksize']].to_csv(output,encoding='utf-8',index=False)

def analysis_frame1(data,output="CIdict_frame01.pickle"):
    temp_sample = ''
    ignore_sample=''
    dictionary={}
    kmers_data = pd.read_csv(data,sep=",")
    for index, row in kmers_data.iterrows():
        if temp_sample!=row["query"]:
            temp_sample = row["query"]
            if row["ksize"] not in dictionary:
                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                if row['frame'] == 'frame_1': 
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
        elif temp_sample==row["query"]:
            if row["ksize"] not in dictionary:
                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                if row['frame'] == 'frame_1': 
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

    with open(output, 'wb') as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

def analysis_frame1(data,output="CIdict_frame01.pickle"):
    temp_sample = ''
    ignore_sample=''
    dictionary={}
    kmers_data = pd.read_csv(data,sep=",")
    for index, row in kmers_data.iterrows():
        if temp_sample!=row["query"]:
            temp_sample = row["query"]
            if row["ksize"] not in dictionary:
                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                if row['frame'] == 'frame_1': 
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
        elif temp_sample==row["query"]:
            if row["ksize"] not in dictionary:
                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
                if row['frame'] == 'frame_1': 
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

    with open(output, 'wb') as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)


def analysis_wrongFrames(data,output="CIdict_wrongframes.pickle"):
    temp_sample = ''
    ignore_sample=''
    dictionary={}
    kmers_data = pd.read_csv(data,sep=",")
    #print(kmers_data.describe())
    correctFrames_df1 = kmers_data[kmers_data["frame"] == "frame_01"][['max_containment','ksize','frame']]
    correctFrames_df1.to_csv('correctFrames_df1.csv')
    wrongFrames_df = kmers_data[kmers_data["frame"] != "frame_01"][['max_containment','ksize','frame']]
    wrongFrames_df.to_csv('wrongFrames_df.csv')
    for index, row in wrongFrames_df.iterrows():
        if row['ksize'] not in dictionary:
            dictionary[row['ksize']] = []
            dictionary[row['ksize']].append(row['max_containment'])
        else:
            dictionary[row['ksize']].append(row['max_containment'])

    print(dictionary)

    with open(output, 'wb') as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

def grab_containment_from_mat(mat_df,ksize):
    """This function converts matrix into df removes pairwise information"""
    """When running contained_by() from sourmash api, a matrix via a csv file called compare.csv is produced"""
    mat = pd.read_csv(mat_df,sep=',')
    df = pd.DataFrame(columns=['A','B','containment','ksize'],dtype=int)
    for i in range(len(mat)):
        for j in range(len(mat)):
            if i > j:
                containment = mat.iloc[j,i]
                df = pd.concat([df, pd.DataFrame.from_records([{'A':mat.columns[i],
                'B':mat.columns[j],
                'containment':containment,
                'ksize':ksize}])],join='inner')
    return(df)
