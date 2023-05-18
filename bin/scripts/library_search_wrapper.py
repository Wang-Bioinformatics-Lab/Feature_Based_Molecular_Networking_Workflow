#!/usr/bin/python


import sys
import getopt
import os
import json
import argparse
import uuid
import pandas as pd
from collections import defaultdict

def chunks(l, n):
    # For item i in a range that is a length of l,
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        yield l[i:i+n]

def search_wrapper(search_param_dict):
    search_files(search_param_dict["spectra_files"], search_param_dict["temp_folder"], search_param_dict["tempresults_folder"], search_param_dict["args"], search_param_dict["params_object"], search_param_dict["library_files"])

def search_files(spectrum_file, library_file, temp_folder, tempresults_folder, path_to_convert, path_to_main_exec,
    min_matched_peaks=6, top_k_results=1, ion_tolerance=0.5, pm_tolerance=20, analog_search=0, max_shift_mass=0.5):


    parameter_filename = os.path.join(temp_folder, str(uuid.uuid4()) + ".params")
    output_parameter_file = open(parameter_filename, "w")

    #Search Criteria
    output_parameter_file.write("MIN_MATCHED_PEAKS_SEARCH={}\n".format(min_matched_peaks))
    output_parameter_file.write("TOP_K_RESULTS={}\n".format(top_k_results))
    output_parameter_file.write("search_peak_tolerance={}\n".format(ion_tolerance))
    output_parameter_file.write("search_parentmass_tolerance={}\n".format(pm_tolerance))
    output_parameter_file.write("ANALOG_SEARCH={}\n".format(0))
    output_parameter_file.write("MAX_SHIFT_MASS={}\n".format(1999))

    #Filtering Criteria
    output_parameter_file.write("FILTER_PRECURSOR_WINDOW={}\n".format(1))
    output_parameter_file.write("MIN_PEAK_INT={}\n".format(0))
    output_parameter_file.write("WINDOW_FILTER={}\n".format(1))
    output_parameter_file.write("FILTER_LIBRARY={}\n".format(1))

    #Scoring Criteria
    output_parameter_file.write("SCORE_THRESHOLD={}\n".format(0.7))

    #Output
    output_parameter_file.write("RESULTS_DIR={}\n".format(tempresults_folder))

    output_parameter_file.write("NODEIDX={}\n".format(0))
    output_parameter_file.write("NODECOUNT={}\n".format(1))

    output_parameter_file.write("EXISTING_LIBRARY_MGF={}\n".format(library_file))

    # Formatting and converting query files
    all_query_spectra_list = []
    
    fileName, fileExtension = os.path.splitext(os.path.basename(spectrum_file))
    output_filename = ""

    if spectrum_file.find("mzXML") != -1 or spectrum_file.find("mzxml") != -1 or spectrum_file.find("mzML") != -1:
        output_filename = os.path.join(temp_folder, fileName + ".pklbin")
        cmd = "%s %s %s" % (path_to_convert, spectrum_file, output_filename)
        print(cmd)
        os.system(cmd)
    else:
        output_filename = os.path.join(temp_folder, os.path.basename(spectrum_file))
        cmd = "cp %s %s" % (spectrum_file, output_filename)
        os.system(cmd)

    #Input
    faked_output_filename = os.path.join(temp_folder, os.path.basename(spectrum_file))
    all_query_spectra_list.append(faked_output_filename)

    output_parameter_file.write("searchspectra=%s\n" % (" ".join(all_query_spectra_list)))
    output_parameter_file.close()

    cmd = "{} ExecSpectralLibrarySearchMolecular {} -ccms_input_spectradir {} -ccms_results_prefix {} -ll 9".format(path_to_main_exec, \
        parameter_filename, temp_folder, tempresults_folder)
    print(cmd)
    os.system(cmd)

    #Removing the spectrum
    try:
        os.remove(output_filename)
    except:
        print("Can't remove", output_filename)


def main():
    parser = argparse.ArgumentParser(description='Running library search parallel')
    parser.add_argument('spectrum_file', help='spectrum_file')
    parser.add_argument('library_file', help='library_file')
    parser.add_argument('result_folder', help='output folder for results')
    parser.add_argument('convert_binary', help='conversion binary')
    parser.add_argument('librarysearch_binary', help='librarysearch_binary')

    args = parser.parse_args()

    temp_folder = "temp"
    try:
        os.mkdir(temp_folder)
    except:
        print("folder error")

    tempresults_folder = "tempresults"
    try:
        os.mkdir(tempresults_folder)
    except:
        print("folder error")

    print(args)

    # performing the search
    search_files(args.spectrum_file, args.library_file, \
        temp_folder, tempresults_folder, \
        args.convert_binary, args.librarysearch_binary)

    # Reformatting the output
    output_results_file = os.path.join(args.result_folder, os.path.basename(args.spectrum_file) + "_" + os.path.basename(args.library_file) + ".tsv")
    
    results_df = pd.read_csv(os.path.join(tempresults_folder, "tempresults"), sep="\t")
    results_df.to_csv(output_results_file, sep="\t", index=False)








if __name__ == "__main__":
    main()
