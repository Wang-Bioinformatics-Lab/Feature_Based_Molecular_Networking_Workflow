#!/usr/bin/python

import sys
import getopt
import os
import json
import argparse
import shutil
import glob
import openms_formatter
import optimus_formatter
import msdial_formatter
import metaboscape_formatter
import xcms_formatter
import mzmine_formatter
import progenesis_formatter
import mztabm_formatter
import agilent_formatter
import agilent2_formatter
from agilent_explorer2_formatter import agilent_explorer2_formatter


def main():
    parser = argparse.ArgumentParser(description='Create parallel parameters')
    parser.add_argument('toolname', help='name of input tool')
    parser.add_argument('quantification_table', help='quantification_table')
    parser.add_argument('quantification_table_reformatted', help='quantification_table_reformatted')
    parser.add_argument('input_spectra_folder', help='input_spectra_folder')
    parser.add_argument('output_mgf', help='output_mgf')
    parser.add_argument('output_metadata', help='Output name of the metadata file, only used for AGILENT_EXPLORER2')
    parser.add_argument('--QUANT_FILE_NORM', default="None", help='QUANT_FILE_NORM')
    
    args = parser.parse_args()

    input_filenames = glob.glob(os.path.join(args.input_spectra_folder, "*"))

    # might be MZMINE MZMINE2 or MZMINE3 in the future
    if "MZMINE" in args.toolname:
        print(args.toolname)

        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        shutil.copyfile(input_mgf, args.output_mgf)
        mzmine_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)

    elif args.toolname == "OPENMS":
        print("OPENMS")
        
        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        shutil.copyfile(input_mgf, args.output_mgf)
        openms_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
    elif args.toolname == "OPTIMUS":
        print("OPTIMUS")

        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        shutil.copyfile(input_mgf, args.output_mgf)
        optimus_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
    elif args.toolname == "MSDIAL":
        print("MSDIAL")
        
        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        shutil.copyfile(input_mgf, args.output_mgf)
        msdial_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
    elif args.toolname == "MSDIAL5":
        print("MSDIAL5")
        
        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        
        import msdial5_formatter
        msdial5_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)

        msdial5_formatter.convert_mgf(input_mgf, args.output_mgf)
    elif args.toolname == "METABOSCAPE":
        print("METABOSCAPE")
        
        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        shutil.copyfile(input_mgf, args.output_mgf)
        metaboscape_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
    elif args.toolname == "XCMS3":
        print("XCMS3")
        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        shutil.copyfile(input_mgf, args.output_mgf)
        xcms_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
    elif args.toolname == "PROGENESIS":
        print("PROGENESIS")

        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]

        compound_scan_mapping = progenesis_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
        progenesis_formatter.convert_mgf(input_mgf, args.output_mgf, compound_scan_mapping)
    elif args.toolname == "AGILENT":
        print("AGILENT")
        
        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        input_mgf = input_filenames[0]
        agilent_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
        agilent_formatter.convert_mgf(input_mgf, args.output_mgf)
    elif args.toolname == "AGILENT2":
        print("AGILENT2")

        feature_table_df = agilent2_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
        agilent2_formatter.convert_mgf(input_filenames, feature_table_df, args.output_mgf)
    elif args.toolname == "AGILENT_EXPLORER2":
        print("AGILENT_EXPLORER2")

        if len(input_filenames) != 1:
            print("Must input exactly 1 spectrum mgf file")
            exit(1)

        # Parse the file
        feature_table_df, metadata, ms2 = agilent_explorer2_formatter.process_pfa_file(input_filenames[0], response_type="Algo", unzipped_files_dir="./PFA_unzipped_files")

        # Write the outputs
        agilent_explorer2_formatter.convert_to_feature_csv(feature_table_df, args.quantification_table_reformatted)
        agilent_explorer2_formatter.write_to_mgf(args.output_mgf, ms2)
        agilent_explorer2_formatter.write_metadata(args.output_metadata, metadata)
    
    
    elif args.toolname == "MZTABM":
        print("MZTABM")

        compound_filename_mapping = mztabm_formatter.convert_to_feature_csv(args.quantification_table, args.quantification_table_reformatted)
        # TODO: Fix This
        mztabm_formatter.create_mgf(input_filenames, args.output_mgf, compound_filename_mapping, name_mangle_mapping=name_mangle_mapping)

    # Finally, we can renormlize the output
    try:
        if args.QUANT_FILE_NORM == "RowSum":
            import pandas as pd
            quant_df = pd.read_csv(args.quantification_table_reformatted, sep=",")
            quant_df = quant_df.loc[:, ~quant_df.columns.str.contains('^Unnamed')]

            for column in quant_df:
                if "Peak area" in column:
                    quant_df[column] = quant_df[column] / sum(quant_df[column]) * 1000000

            quant_df.to_csv(args.quantification_table_reformatted, sep=",", index=False)
    except:
        pass

if __name__ == "__main__":
    main()
