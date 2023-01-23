import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import warnings

def CIbox_frame01(data="dictionary.pickle",wd="data/"):
    with open(data, 'rb') as handle:
        unserialized_data = pickle.load(handle)
    warnings.filterwarnings('ignore')
    #box plot #1

    df1 = pd.DataFrame.from_dict(unserialized_data).T.explode('highest')
    df1.reset_index(inplace=True)
    df1 = df1.set_index([df1.index, 'index'])['highest'].unstack()
    df1 = df1.apply(lambda x: pd.Series(x.dropna().values))

    vals1, names1, xs1 = [],[],[]
    for i, col in enumerate(df1.columns):
        vals1.append(df1[col].values)
        names1.append(str(col)+"-mer")
        xs1.append(np.random.normal(i + 1, 0.04, df1[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

    # box plot #2

    df2 = pd.DataFrame.from_dict(unserialized_data).T.explode('2ndhighest')
    df2.reset_index(inplace=True)
    df2 = df2.set_index([df2.index, 'index'])['2ndhighest'].unstack()
    df2 = df2.apply(lambda x: pd.Series(x.dropna().values))

    vals2, names2, xs2 = [],[],[]
    for i, col in enumerate(df2.columns):
        vals2.append(df2[col].values)
        names2.append(str(col)+"-mer")
        xs2.append(np.random.normal(i + 1, 0.04, df2[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))
    bplot1 = ax1.boxplot(vals1, labels=names1, notch=False, showmeans=False)
    bplot2 = ax2.boxplot(vals2, labels=names2, notch=False, showmeans=False)

    palette = ['red', 'gray', 'blue', 'yellow','orange','purple','brown','pink','olive','cyan']

    for xA, xB, valA, valB, c in zip(xs1, xs2, vals1, vals2, palette):
        ax1.scatter(xA, valA, alpha=0.4, color=c)
        ax2.scatter(xB, valB, alpha=0.4, color=c)

    sns.set_style("whitegrid")
    ind = np.arange(11) 
    for ax in fig.get_axes():
        ax.grid(visible=None)
        ax.set_xticks(ticks=ind,labels=['','k7','k14','k21','k28','k35','k42','k49','k56','k63','k70'],fontsize=15)
        ax.set_ylim(0,1.1)

    fig.suptitle("Highest andf Second Highest Containment Indexes\nReported for the Correct ORFs Across Kmers by Sourmash",fontsize=20)
    fig.text(0.07,0.5,'Containment Index',ha='center',va='center',rotation='vertical',fontsize=15)

    fig.savefig(wd+"boxplot.jpeg",bbox_inches='tight')

def CIbox_wrongORFs(data="dictionary.pickle",wd="data/"):
    with open(data, 'rb') as handle:
        unserialized_data = pickle.load(handle)
    warnings.filterwarnings('ignore')

    df1 =pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in unserialized_data.items() ]))
    df1 = pd.DataFrame.from_dict(df1)
    plt.boxplot(df1)
    plt.savefig("test_wrong.jpeg")
#    df1.reset_index(inplace=True)
#    print(df1.head())
    #df1.reset_index(inplace=True)
    #df1 = df1.set_index([df1.index, 'index'])['highest'].unstack()
    #df1 = df1.apply(lambda x: pd.Series(x.dropna().values))

def CIhist(input="containment.csv",kmers=[7,14,21,28,35,42,49,56,63,70],wd="results/"):
    data = pd.read_csv(input,sep=",")

    def frame(x):
        return(x[-7:])

    match_name=1
    max_containment=0
    labels=[]
    frame_1_lst=[]
    frame_x_lst=[]

    for i in kmers:
        labels.append(str(i)+"-mer")
        df = data[data["ksize"] == i]
        df["ksize"] = pd.to_numeric(df["ksize"])
        df = df.sort_values(by=['query', "ksize", "max_containment"],ascending=False)
        df = df.groupby(["query","ksize"]).head(1).reset_index(drop=True)
        frame_1_containment = df['frame'].value_counts()[0]

        if i < 21: #kmers 21 and above report only containment index of 1
            frame_x_containment = df['frame'].value_counts()[1]
        else:
            frame_x_containment = 0
        total = frame_1_containment+frame_x_containment
        frame_1_lst.append(frame_1_containment/total*100)
        frame_x_lst.append(frame_x_containment/total*100)

        # make data:
        x = np.array(["Correct", "Incorrect"])
        y = np.array([frame_1_containment/total*100,frame_x_containment/total*100])


    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots()
    #plt.figure(figsize=(20,10))

    rects1 = ax.bar(x - width/2, frame_1_lst, width, label='Correct', color="darkturquoise",edgecolor='darkturquoise')
    rects2 = ax.bar(x + width/2, frame_x_lst, width, label='Incorrect', color="greenyellow",edgecolor="greenyellow")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('%',fontsize = 15)
    ax.set_title('Percentage of Correct ORF Report by K-size',fontsize = 20)

    ax.set_xticks(x)
    ax.set_xticklabels(labels,fontsize=15,rotation=45)
    ax.set_yticklabels(labels=[0.0,0.2,0.4,0.6,0.8,1.0],fontsize=12)

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize=15)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.3)
    ax.spines['left'].set_linewidth(0.3)

    ax.grid(visible=None)

    fig.tight_layout()
    fig.set_size_inches(18.5, 10.5)
    plt.savefig(wd+'Percentage_Correct_ORF_K-size.jpeg',bbox_inches='tight')
