import pandas as pd
import csv
import os, glob, subprocess, argparse
from frame_predict import report_containment

def main(args):
    #arguments
    wd = args.wd #working directory where output will go
    ksizes = [int(item) for item in args.ksizes.split(',')] # returns list of ksizes

    #Getting containment index that we are interested in. Containment index obtained from comparing all queries to the reference
    for ksize in ksizes:
        """Extract containment indexes from matrix for each ksize sourmash comapre result"""
        df=report_containment.prep_df_for_compare_containment_Ksize(compare_K_csv_file=wd+'compare'+str(ksize)+'.csv',ksize=ksize) #extract containment indexes from matrix for each ksize sourmash comapre result
        df.to_csv(wd+'compare_containment'+str(ksize)+'.csv',sep=',') 
        df_2nd_highest=report_containment.report_df_of_second_highest_containment_indexes(ksize=ksize,data=wd+'compare_containment'+str(ksize)+'.csv')
        df_2nd_highest.to_csv(wd+'second_highest_containment'+str(ksize)+'.csv')

    #output all csv containment index files together by forst removing all headers except from the first file 
     #Produce output file to contain header (once)
    #FNR represents the number of the processed record in a single file. And NR represents it globally, so first line is accepted and the rest are ignored as before
    cmd = f"awk '(NR == 1) || (FNR > 1)' {wd}compare_containment*.csv | grep -v 'uniprotkb.fasta' > {wd}containment.csv"
    subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)

    #output all csv containment index files together by forst removing all headers except from the first file for the second highest containment indexes report
    cmd1 = f"awk '(NR == 1) || (FNR > 1)' {wd}second_highest_containment*.csv > {wd}second_highest_containment.csv"
    subprocess.run(cmd1, stdout=subprocess.PIPE, shell=True)

    #### report CI for frames 1
    frame_1_csv=report_containment.report_df_of_frame1(wd+"containment.csv")
    frame_1_csv.to_csv(wd+"frame1_containment.csv",encoding='utf-8',index=False)

    #### report CI for all other frames
    frame_X_csv=report_containment.report_df_of_frameX(wd+"containment.csv")
    frame_X_csv.to_csv(wd+'frameX_containment.csv',encoding='utf-8',index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Anayze whether using minhash containment index via Sourmash propgram can accurately predict the correct reading frame') 

    parser.add_argument('--ksizes', default='7,14,21', help='K-mer size for containment index, a parameter used in Sourmash sketch.', type=str)

    parser.add_argument('--wd', help = 'Output files from sourmash program')

    args = parser.parse_args()

    main(args)

