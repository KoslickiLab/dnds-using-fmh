import pickle, os, glob, subprocess
import pandas as pd
import sourmash as sm

#Creating signature of my ref sig to be the same scaled value of my query sig
#query_input='/data/jzr5814/data/frame_analysis_using_bash_script/query_frames.faa'
#scaled=100
#query_output='/data/jzr5814/data/frame_analysis_using_bash_script/query_frames_scaled100.sig.zip'
#cmd2=f"sourmash sketch protein {query_input} -p k=7,k=14,k=21,k=28,k=35,k=42,k=49,k=56,k=63,k=70,scaled={scaled} --singleton -o {query_output}"
#print(cmd2)
#subprocess.run(cmd2, stdout=subprocess.PIPE, shell=True)

#output containment csv file
query_file = '/data/jzr5814/data/frame_analysis_using_bash_script/query_frames_scaled100.sig.zip'
ref_file = '/data/jzr5814/data/uniprotkb.sig'

ksize_lst = [7,14,21,28,35,42,49,56,63,70]
 
df=pd.DataFrame(columns=['query','frame','ksize','containment'])

for k in ksize_lst:
    query_sigs = sm.load_file_as_signatures(query_file, ksize=k)
    ref_sigs = sm.load_file_as_signatures(ref_file, ksize=k)

    for sig in query_sigs:
        for sig2 in ref_sigs:
            c = sig.contained_by(sig2)
                    
        df=df.append({'query':sig,'frame':str(sig)[-1:],'ksize':k,'containment':c},ignore_index=True)               

df.to_csv('/data/jzr5814/data/frame_analysis_using_bash_script/containment.csv')
