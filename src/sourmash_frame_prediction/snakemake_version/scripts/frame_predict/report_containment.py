import pickle, os, glob, subprocess
import pandas as pd

#def produce_containment_csv(wd):
#    cmd1 = f"awk '(NR == 1) || (FNR > 1)' {wd}prefetch_res_*.csv | cut -d ',' -f 3,7,14,15 | tr '_' ',' | sed 's/max,containment/max_containment/' > {wd}results_kmers.csv"
#    subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

def extract_frame_number(x):
    """Extract the frame number from the fasta sequence name
    When creating the six reading frame .faa file,
    the fasta names are modified to remove special characters 
    and to indicate the frame number that the sequence represents by adding '_X' (X being a frame number)
    This function takes in that fasta sequence name and returns frame_X."""
    return('frame_'+x[-1:])

#def prep(data,output="containment.csv"):
#    kmers_data = pd.read_csv(data,sep=",",converters={"name":str}).rename(columns = {'match':'query','name':'frame'}).sort_values(by=["query", "ksize", "max_containment"],ascending=False)
    #match_name=1
    #max_containment=0
#    kmers_data["frame"] = kmers_data["frame"].astype(str).apply(frame)

 #   kmers_data.to_csv(output,encoding='utf-8')

def prep_df_for_compare_containment_Ksize(compare_K_csv_file, ksize):
    """Return datafrane of prepared compare_containment_K.csv file.
    This function takes in the csv output from sourmash compare and 
    prepares it for further analyses (figures). Because this csv file is a matrix,
    this function transposes the matrix to keep the column names (to indicare the query being compared to the reference) 
    and the first row (containment values of queries to the reference). A frame column that indicates frame number is added"""
    df = pd.read_csv(compare_K_csv_file, sep=',')
    transpose_df=df.T.iloc[:,0].reset_index().rename(columns={'index':'query',0:'containment'}).rename_axis('index') #extract column where all queries are compared to ref
    transpose_df['ksize'] = ksize
    transpose_df["frame"] = transpose_df["query"].astype(str).apply(extract_frame_number)
    return(transpose_df)

def report_df_of_frame1(data="containment.csv",output="frame1_containment.csv"):
    """Return csv file that only reports containment indexes that involve sourmash comapre of frame 1, which is the correct fraame."""
    frame1_data = pd.read_csv(data,sep=",").reset_index()
    frame1_data = frame1_data[frame1_data['frame'] == 'frame_1']
    return(frame1_data[['index','containment','frame','ksize']])
#    frame1_data[['index','containment','frame','ksize']].to_csv(output,encoding='utf-8',index=False)
#    frame1_data[['index','max_containment','frame','ksize']].to_csv(output,encoding='utf-8',index=False) #columns based off of prefetch results

def report_df_of_frameX(data="containment.csv",output="frameX_containment.csv"):
    """Return csv file that only reports containment indexes that involve sourmash comapre of all frames except for frame 1, which are the incorrect frames."""
    framex_data = pd.read_csv(data,sep=",").reset_index()
    framex_data = framex_data[framex_data['frame'] != 'frame_1']
    return(framex_data[['index','containment','frame','ksize']])
#    framex_data[['index','containment','frame','ksize']].to_csv(output,encoding='utf-8',index=False)
#    framex_data[['index','max_containment','frame','ksize']].to_csv(output,encoding='utf-8',index=False) #columns based off of prefetch results

def report_df_of_second_highest_containment_indexes(data=compare_containment_K_csv_file,ksize):
    """"""
    df=pd.read_csv(compare_containment_K_csv_file,sep=',')[['query','ksize','containment','frame']].iloc[1:]
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
    return(pd.DataFrame.from_dict(new_dict).rename_axis('index'))

#def analysis_frame1(data,output="CIdict_frame01.pickle"):
#    temp_sample = ''
#    ignore_sample=''
#    dictionary={}
#    kmers_data = pd.read_csv(data,sep=",")
#    for index, row in kmers_data.iterrows():
#        if temp_sample!=row["query"]:
#            temp_sample = row["query"]
#            if row["ksize"] not in dictionary:
#                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    if ignore_sample != temp_sample:
#                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#                        ignore_sample=temp_sample
#            elif row["ksize"] in dictionary:
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    if ignore_sample != temp_sample:
#                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#                        ignore_sample=temp_sample
#        elif temp_sample==row["query"]:
#            if row["ksize"] not in dictionary:
#                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#            elif row["ksize"] in dictionary:
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    if ignore_sample != temp_sample:
#                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#                        ignore_sample=temp_sample

#    with open(output, 'wb') as handle:
#        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)

#def analysis_frame1(data,output="CIdict_frame01.pickle"):
#    temp_sample = ''
#    ignore_sample=''
#    dictionary={}
#    kmers_data = pd.read_csv(data,sep=",")
#    for index, row in kmers_data.iterrows():
#        if temp_sample!=row["query"]:
#            temp_sample = row["query"]
#            if row["ksize"] not in dictionary:
#                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    if ignore_sample != temp_sample:
#                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#                        ignore_sample=temp_sample
#            elif row["ksize"] in dictionary:
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    if ignore_sample != temp_sample:
#                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#                        ignore_sample=temp_sample
#        elif temp_sample==row["query"]:
#            if row["ksize"] not in dictionary:
#                dictionary[row['ksize']] = {'highest':[],'2ndhighest':[]}
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#            elif row["ksize"] in dictionary:
#                if row['frame'] == 'frame_1': 
#                    dictionary[row['ksize']]['highest'].append(row['max_containment'])
#                else:
#                    if ignore_sample != temp_sample:
#                        dictionary[row['ksize']]['2ndhighest'].append(row['max_containment'])
#                        ignore_sample=temp_sample

#    with open(output, 'wb') as handle:
#        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)


#def analysis_wrongFrames(data,output="CIdict_wrongframes.pickle"):
#    temp_sample = ''
#    ignore_sample=''
#    dictionary={}
#    kmers_data = pd.read_csv(data,sep=",")
    #print(kmers_data.describe())
#    correctFrames_df1 = kmers_data[kmers_data["frame"] == "frame_01"][['max_containment','ksize','frame']]
#    correctFrames_df1.to_csv('correctFrames_df1.csv')
#    wrongFrames_df = kmers_data[kmers_data["frame"] != "frame_01"][['max_containment','ksize','frame']]
#    wrongFrames_df.to_csv('wrongFrames_df.csv')
#    for index, row in wrongFrames_df.iterrows():
#        if row['ksize'] not in dictionary:
#            dictionary[row['ksize']] = []
#            dictionary[row['ksize']].append(row['max_containment'])
#        else:
#            dictionary[row['ksize']].append(row['max_containment'])

#    print(dictionary)

#    with open(output, 'wb') as handle:
#        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
