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

    # These are params to be able to merge things appropriately
    parser.add_argument('--topk', default=None, help='When merging, we want to include the top k results for each unique key')
    parser.add_argument('--key_column', default=None, help='column of the key to group by')
    parser.add_argument('--sort_column', default=None, help='column of the to sort by')

    args = parser.parse_args()

    all_results_files = glob.glob(os.path.join(args.input_folder, "*.tsv"))

    all_results_list = []
    for i, results_file in enumerate(all_results_files):
        temp_df = pd.read_csv(results_file, sep="\t")

        if len(temp_df) > 0:
            all_results_list.append(temp_df)
    
    # merging results
    all_results_df = pd.concat(all_results_list, ignore_index=True)

    # Filtering when appropriate
    if args.topk is not None:
        topk_filter = int(args.topk)

        all_results_df = all_results_df.sort_values(by=args.sort_column, ascending=False)
        all_results_df = all_results_df.groupby(args.key_column).head(int(args.topk))

    # writing results
    all_results_df.to_csv(args.output_file, sep="\t", index=False)

if __name__ == "__main__":
    main()
