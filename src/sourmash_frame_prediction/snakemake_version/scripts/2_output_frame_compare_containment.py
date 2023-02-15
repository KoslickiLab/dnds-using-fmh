import pandas as pd
import csv
import os, glob, subprocess, argparse
from frame_predict import report_containment

def main(args):
    #wd = '/data/jzr5814/data/frame_analysis_using_bash_script/compare_by_multiple_jobs/test1/'
    #ksizes = [7,14,21,28,35,42,49,56,63,70]
    wd = args.wd #working directory where output will go
    ksizes = [int(item) for item in args.ksizes.split(',')] # returns list of ksizes

    #Getting containment index that we are interested in. Containment index obtained from comparing all queries to the reference
    for ksize in ksizes:
        """Extract containment indexes from matrix for each ksize sourmash comapre result"""
        df=prep_df_for_compare_containment_Ksize(wd+'compare'+str(ksize)+'.csv',ksize) #extract containment indexes from matrix for each ksize sourmash comapre result
        df.to_csv(wd+'compare_containment'+str(ksize)+'.csv',sep=',') 
        df_2nd_highest=report_df_of_second_highest_containment_indexes(wd+'compare_containment'+str(ksize)+'.csv')
        df_2nd_highest.to_csv(wd+'second_highest_containment'+str(ksize)+'.csv')

    cmd = f"awk '(NR == 1) || (FNR > 1)' {wd}compare_containment*.csv | grep -v 'uniprotkb.fasta' > {wd}containment.csv"
    subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)

    cmd1 = f"awk '(NR == 1) || (FNR > 1)' {wd}second_highest_containment*.csv > {wd}second_highest_containment.csv"
    subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

    #### report CI for frames 1
    frame_1_csv=extract_frame1(wd+"containment.csv",{wd}+"frame1_containment.csv")
    frame_1_csv.to_csv(output,encoding='utf-8',index=False)

    #### report CI for all other frames
    frame_X_csv=extract_frameX(wd+"containment.csv",{wd}+"frameX_containment.csv")
    frame_X_csv.to_csv('frameX_containment.csv',encoding='utf-8',index=False)

    #create figures of CI analysis
    figures.CIbox_frames(frame_1data=wd+"frame1_CI.csv",frame_xdata=wd+"framex_CI.csv",output=wd+'boxplot_frame1_and_all_other_frames.jpeg') #all other frames not 1
    figures.CIbox_frames(frame_1data=wd+"frame1_CI.csv",frame_xdata=wd+"second_highest_containment.csv",output=wd+'boxplot_frame1_and_second_highest.jpeg') #second highes

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Anayze whether using minhash containment index via Sourmash propgram can accurately predict the correct reading frame') 

    parser.add_argument('--ksizes', default='7,14,21', help='K-mer size for containment index, a parameter used in Sourmash sketch.', type=str)

    parser.add_argument('--wd', help = 'Output files from sourmash program')

    args = parser.parse_args()

    main(args)

