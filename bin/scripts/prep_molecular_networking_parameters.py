#!/usr/bin/python


import sys
import getopt
import os
import json
import argparse

def number_scans_in_mgf_file(mgf_filename):
    number_spectra = 0
    number_real_spectra = 0
    number_of_peaks = 0
    for line in open(mgf_filename):
        if line.rstrip() == "BEGIN IONS":
            number_of_peaks = 0
            number_spectra += 1
        if len(line) > 2:
            if line[0].isdigit():
                number_of_peaks += 1
        if line.rstrip() == "END IONS":
            if number_of_peaks > 0:
                number_real_spectra += 1

    return number_spectra, number_real_spectra

def main():
    parser = argparse.ArgumentParser(description='Create parallel parameters')
    parser.add_argument('mgf_filename', help='Input mgf file to network')
    parser.add_argument('parameters_output_folder', help='output folder for parameters')
    parser.add_argument('--parallelism', default=1, type=int, help='Parallelism')
    parser.add_argument('--min_matched_peaks', default=6, type=int, help='min_matched_peaks')
    parser.add_argument('--ms2_tolerance', default=0.5, type=float, help='ms2_tolerance')
    parser.add_argument('--pm_tolerance', default=0.5, type=float, help='pm_tolerance')
    parser.add_argument('--min_cosine', default=0.5, type=float, help='min_cosine')
    parser.add_argument('--max_shift', default=200, type=float, help='max_shift')
    args = parser.parse_args()

    #Determing number of spectra in mgf file
    number_of_spectra, number_real_spectra = number_scans_in_mgf_file(args.mgf_filename)

    parallelism = args.parallelism
    if parallelism > number_of_spectra:
        parallelism = 1

    recommended_parallelism = max(1, int(number_real_spectra/1000))

    print("recommended_parallelism", recommended_parallelism)

    parallelism = min(recommended_parallelism, parallelism)

    number_per_partition = int(number_of_spectra/parallelism)
    for i in range(parallelism):
        output_parameter_file = open(os.path.join(args.parameters_output_folder, str(i) + ".params"), "w")
        output_parameter_file.write("ALIGNS_FORMAT={}\n".format("tsv"))
        output_parameter_file.write("MIN_MATCHED_PEAKS={}\n".format(args.min_matched_peaks))
        output_parameter_file.write("TOLERANCE_PEAK={}\n".format(args.ms2_tolerance))
        output_parameter_file.write("TOLERANCE_PM={}\n".format(args.pm_tolerance))
        output_parameter_file.write("PAIRS_MIN_COSINE={}\n".format(args.min_cosine))
        output_parameter_file.write("MAX_SHIFT={}\n".format(args.max_shift))
        output_parameter_file.write("INPUT_SPECTRA_MS2={}\n".format(args.mgf_filename))

        start_idx = number_per_partition * i
        end_idx = number_per_partition * (i + 1) - 1
        if i == parallelism - 1:
            end_idx = number_of_spectra

        output_parameter_file.write("IDX_START=%d\n" % (start_idx))
        output_parameter_file.write("IDX_END=%d\n" % (end_idx))


if __name__ == "__main__":
    main()
