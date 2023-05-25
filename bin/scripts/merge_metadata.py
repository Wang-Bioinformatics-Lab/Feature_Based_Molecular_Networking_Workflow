import os
import pandas as pd
import argparse

def _load_metadata(input_filename):
    # look at extension
    if input_filename.endswith(".tsv"):
        input_df = pd.read_csv(input_filename, sep="\t")
    elif input_filename.endswith(".csv"):
        input_df = pd.read_csv(input_filename, sep=",")
    elif input_filename.endswith(".xlsx"):
        input_df = pd.read_excel(input_filename)

    return input_df

def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input_metadata')
    parser.add_argument('merged_metadata')
    args = parser.parse_args()

    if not os.path.exists(args.input_metadata):
        # This is likely not valid, lets skip it
        input_metadata = pd.DataFrame()
    else:
        input_metadata = _load_metadata(args.input_metadata)

    # outputing metadata
    input_metadata.to_csv(args.merged_metadata, sep="\t", index=False)


if __name__ == '__main__':
    main()