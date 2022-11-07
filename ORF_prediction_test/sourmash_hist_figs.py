import pickle
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

with open('../data/dictionary.pickle', 'rb') as handle:
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

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.title("Containment Index of Correct Open Reading Frame Reported by Sourmash",fontsize=20)
plt.ylabel('Containment Index',fontsize=15)

plt.savefig('../data/Containment_Index_Correct_ORF.jpeg',bbox_inches='tight')

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

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.title("Second Highest Containment Index Reported by Sourmash",fontsize=20)
plt.ylabel('Containment Index',fontsize=15)

plt.savefig('../data/Second_Highest_Containment_Index_Reported_SourmashF.jpeg',bbox_inches='tight')


