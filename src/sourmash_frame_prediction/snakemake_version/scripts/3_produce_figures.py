import pandas as pd
import csv
import os, glob, subprocess, argparse
from frame_analysis import figures

def main(args):

    frame1_input=args.frame_1data
    frame2_input=args.frame_Xdata
    output = args.output #working directory where output will go

    #create figures of CI analysis
    figures.CIbox_frames(frame_1data=frame1_input,frame_xdata=frame2_input,output=output) #all other frames not 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Produce figures to anayze whether using minhash containment index via Sourmash propgram can accurately predict the correct reading frame') 

    parser.add_argument('--frame_1data', default='frame1_CI.csv', help='', type=str)

    parser.add_argument('--frame_Xdata', default='framex_CI.csv', help='', type=str)

    parser.add_argument('--output', default='boxplot.jpeg', help='Output figure for containment analysis of frames', type=str)

    args = parser.parse_args()

    main(args)

