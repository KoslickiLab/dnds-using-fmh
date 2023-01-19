import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def CIhist(data="dictionary.pickle",wd="data/"):
    with open(data, 'rb') as handle:
        unserialized_data = pickle.load(handle)
        
    df = pd.DataFrame.from_dict(unserialized_data).T.explode('highest')
    df.reset_index(inplace=True)
    df = df.set_index([df.index, 'index'])['highest'].unstack()
    df = df.apply(lambda x: pd.Series(x.dropna().values))

    vals, names, xs = [],[],[]
    for i, col in enumerate(df.columns):
        vals.append(df[col].values)
        names.append(str(col)+"-mer")
        xs.append(np.random.normal(i + 1, 0.04, df[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

    plt.figure(figsize=(18,10))

    plt.boxplot(vals, labels=names)

    palette = ['red', 'gray', 'blue', 'yellow','orange','purple','brown','pink','olive','cyan']

    for x, val, c in zip(xs, vals, palette):
        plt.scatter(x, val, alpha=0.4, color=c)

    plt.grid(visible=None)
    plt.xticks(rotation=45,fontsize=15)
    plt.yticks(fontsize=12)
    plt.ylim(0,1.1)

    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)

    plt.title("Highest Containment Index Across Kmers Reported by Sourmash",fontsize=20)
    plt.ylabel('Containment Index',fontsize=15)

    plt.savefig(wd+'Highest_Containment_Index_Correct_ORF.jpeg',bbox_inches='tight')

    df = pd.DataFrame.from_dict(unserialized_data).T.explode('2ndhighest')
    df.reset_index(inplace=True)
    df = df.set_index([df.index, 'index'])['2ndhighest'].unstack()
    df = df.apply(lambda x: pd.Series(x.dropna().values))

    vals, names, xs = [],[],[]
    for i, col in enumerate(df.columns):
        vals.append(df[col].values)
        names.append(str(col)+"-mer")
        xs.append(np.random.normal(i + 1, 0.04, df[col].values.shape[0]))  # adds jitter to the data points - can be adjusted

    plt.figure(figsize=(18,10))

    plt.boxplot(vals, labels=names)

    palette = ['red', 'gray', 'blue', 'yellow','orange','purple','brown','pink','olive','cyan']

    for x, val, c in zip(xs, vals, palette):
        plt.scatter(x, val, alpha=0.4, color=c)

    plt.grid(visible=None)
    plt.xticks(rotation=45,fontsize=15)
    plt.yticks(fontsize=12)
    plt.ylim(0,1.1)

    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['top'].set_visible(False)

    plt.title("Second Highest Containment Index Across Kmers Reported by Sourmash",fontsize=20)
    plt.ylabel('Containment Index',fontsize=15)

    plt.savefig(wd+"Second_Highest_Containment_Index_Reported_Sourmash.jpeg",bbox_inches='tight')

def CIbox(input="containment.csv",kmers=[7,14,21,28,35,42,49,56,63,70],wd="results/"):
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
