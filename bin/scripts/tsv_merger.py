#!/usr/bin/python


import sys
import getopt
import os
import pandas as pd
from collections import defaultdict
import argparse
import glob


def main():
    # Parsing the arguments
    parser = argparse.ArgumentParser(description='Merging Results Files')
    parser.add_argument('input_folder', help='input_folder')
    parser.add_argument('output_file', help='output_file')

    args = parser.parse_args()

    all_results_files = glob.glob(os.path.join(args.input_folder, "*.tsv"))

    all_results_list = []
    for i, results_file in enumerate(all_results_files):
        temp_df = pd.read_csv(results_file, sep="\t")

        if len(temp_df) > 0:
            all_results_list.append(temp_df)
    
    # merging results
    all_results_df = pd.concat(all_results_list, ignore_index=True)

    # writing results
    all_results_df.to_csv(args.output_file, sep="\t", index=False)

if __name__ == "__main__":
    main()
