#!/usr/bin/env python3
"""Create arguments to run on multiple ksizes for frame analysis"""

import argparse, pickle, os, glob, subprocess
import pandas as pd
import sourmash as sm

def main(args):

    k = args.ksize

    query_file = '/data/jzr5814/data/frame_analysis_using_bash_script/query_frames_scaled100.sig.zip'
    print('load'+query_file)

    ref_file = '/data/jzr5814/data/uniprotkb.sig'
    print('load'+ref_file)

    print('set klist') 
    print(k)

    print('empty df')
    df=pd.DataFrame(columns=['query','frame','ksize','containment'])

    query_sigs = sm.load_file_as_signatures(query_file, ksize=k)
    ref_sigs = sm.load_file_as_signatures(ref_file, ksize=k)

    for sig in query_sigs:
        for sig2 in ref_sigs:
            c = sig.contained_by(sig2)

        print(sig,str(sig)[-1:],k,c)
        df=pd.DataFrame.from_records([{'query':sig,'frame':str(sig)[-1:],'ksize':k,'containment':c}])

    print('saving containment report to csv file')
    df.to_csv('/data/jzr5814/data/frame_analysis_using_bash_script/contained_by_multiple_jobs/containment'+str(k)+'.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = 'obtain containment index between ref and query signatures'
    )

    parser.add_argument(
        '--ksize',
        type=int,
        help = 'ksize'
    )

    args = parser.parse_args()

    main(args)
