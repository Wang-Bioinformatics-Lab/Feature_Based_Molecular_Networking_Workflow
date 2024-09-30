#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mingxun wang
@purpose: to take the quant file and raw data and create a tall table for visualization and linkout
"""
import pandas as pd
import os
import sys
import argparse
import glob

def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='quant file and raw data and create a tall table for visualization and linkout')

    parser.add_argument('quant_file', type=str, help='quant file')
    parser.add_argument('raw_data', type=str, help='raw data folder')

    # output
    parser.add_argument('output_tall', type=str, help='output file')

    args = parser.parse_args()

    # read in the quant file
    quant_df = pd.read_csv(args.quant_file, sep=',')

    # get the raw files
    raw_files = glob.glob(os.path.join(args.raw_data, '*'))

    # getting the filename headers
    all_filenames = []
    for header in quant_df.columns:
        if "Peak area" in header:
            filename = header.replace("Peak area", "").rstrip()
            all_filenames.append(filename)
    

    # Cross check that these are in the raw filesnames
    intesection_filenames = {}
    
    for raw_file in raw_files:
        for filename in all_filenames:
            if os.path.basename(raw_file) == filename:
                intesection_filenames[filename] = raw_file
                break
    
    print(len(intesection_filenames), "raw files found")


    # creating the output tall file, using pandas melt
    all_filenames_columns = ["{} Peak area".format(filename) for filename in all_filenames]
    tall_df = quant_df.melt(id_vars=['row ID'], value_vars=all_filenames_columns, var_name='filename', value_name='area')

    # clean the filename
    tall_df['filename'] = tall_df['filename'].str.replace(" Peak area", "")
    
    # intersect with a 0 or 1, if the raw data is present
    tall_df['raw_data'] = 0
    for filename in intesection_filenames.keys():
        tall_df.loc[tall_df['filename'] == filename, 'raw_data'] = 1

    # Adding in the m/z and retention time by merging from the original quant_df
    tall_df = pd.merge(tall_df, quant_df[['row ID', 'row m/z', 'row retention time']], on='row ID', how='left')

    # Write out
    tall_df.to_csv(args.output_tall, index=False, sep="\t")


if __name__=="__main__":
    # there should be obly one input file
    main()
