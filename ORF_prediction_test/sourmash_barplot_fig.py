import pickle
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

data = pd.read_csv("../data/results_kmers.csv",sep=",",header=None,engine="python",on_bad_lines="skip",quoting=3)

def frame(x):
    return(x[-7:])
def sample(x):
    return(x[:15])

match_name=1
max_containment=0
ksize=2
labels=[]
frame_1_lst=[]
frame_x_lst=[]

for i in [7,14,21,28,35,42,49,56,63,77]:
    labels.append(str(i)+"-mer")
    smallk = data[data[ksize] == i]
    smallk["frame"] = smallk[match_name].apply(frame)
    smallk = smallk[smallk['frame'].str.startswith('frame_')]
    smallk['sample'] = smallk[match_name].apply(sample)
    smallk[ksize] = pd.to_numeric(smallk[ksize])
    smallk = smallk[["frame","sample",max_containment,ksize]].rename(columns={0: 'max_containment', 2: 'ksize'})
    smallk = smallk.sort_values(by=['sample', "ksize", "max_containment"],ascending=False)
    smallk = smallk.groupby(["sample","ksize"]).head(1).reset_index(drop=True)
    frame_1_containment = smallk['frame'].value_counts()[0]
    if i < 21:
        frame_x_containment = smallk['frame'].value_counts()[1]
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
plt.savefig('../data/Percentage_Correct_ORF_K-size.jpeg',bbox_inches='tight')

