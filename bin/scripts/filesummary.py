#!/usr/bin/python


import sys
import getopt
import os
import json
import argparse
import uuid
from collections import defaultdict
#import ming_fileio_library
import ming_parallel_library
import csv
import re
import pandas as pd
import glob

def summary_wrapper(search_param_dict):
    summary_files(search_param_dict["spectrum_file"], search_param_dict["tempresults_folder"], search_param_dict["args"])

def summary_files(spectrum_file, tempresults_folder, args):
    summary_filename = os.path.join(tempresults_folder, os.path.basename(spectrum_file) + ".summary")
    cmd = "export LC_ALL=C && {} {} -x \"run_summary delimiter=tab\" > {}".format(args.msaccess_binary, spectrum_file, summary_filename)
    print(cmd)
    os.system(cmd)

def main():
    parser = argparse.ArgumentParser(description='Running library search parallel')
    parser.add_argument('spectra_folder', help='spectrafolder')
    parser.add_argument('result_file', help='output folder for parameters')
    parser.add_argument('msaccess_binary', help='output folder for parameters')
    parser.add_argument('--parallelism', default=1, type=int, help='Parallelism')
    args = parser.parse_args()

    spectra_files = glob.glob(os.path.join(args.spectra_folder, "*"))
    spectra_files.sort()

    tempresults_folder = "tempresults"
    try:
        os.mkdir(tempresults_folder)
    except:
        print("folder error")

    parameter_list = []
    for spectrum_file in spectra_files:
        param_dict = {}
        param_dict["spectrum_file"] = spectrum_file
        param_dict["tempresults_folder"] = tempresults_folder
        param_dict["args"] = args

        parameter_list.append(param_dict)

    print("Parallel to execute", len(parameter_list))
    ming_parallel_library.run_parallel_job(summary_wrapper, parameter_list, 10)


    """Merging Files and adding full path"""
    all_result_files = glob.glob(os.path.join(tempresults_folder, "*"))
    full_result_list = []
    for input_file in all_result_files:
        try:
            results_df = pd.read_csv(input_file, sep="\t")
            result_list = results_df.to_dict(orient="records")
            for result in result_list:
                output_dict = {}
                output_dict["Filename"] = result["Filename"]
                output_dict["Vendor"] = result["Vendor"]
                output_dict["Model"] = result["Model"]
                output_dict["MS1s"] = result["MS1s"]
                output_dict["MS2s"] = result["MS2s"]
                full_result_list.append(output_dict)
        except:
            #raise
            print("Error", input_file)

    pd.DataFrame(full_result_list).to_csv(args.result_file, sep="\t", index=False)



if __name__ == "__main__":
    main()
