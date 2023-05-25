#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mingxun wang
@purpose: to convert the Agilent file into a diserable format
"""
import pandas as pd
import sys
from pyteomics import mgf
from collections import defaultdict
import ming_spectrum_library


def _convert_details_feature_csv(input_filename, output_filename):
    input_df = pd.read_csv(input_filename, sep=',', skiprows=2, index_col=False)

    # We need to turn this table into a proper table
    input_df = input_df[["Label", "Mass", "RT", "Area", "File"]]
    input_df["row ID"] = input_df["Label"].apply(lambda x: x.split(" ")[1].split(":")[0])

    # Dropping the extension of the filename and keeping only the basename
    input_df["File"] = input_df["File"].apply(lambda x: x.split("/")[-1].split(".")[0])

    # Adding Peak area suffix
    input_df["File"] = input_df["File"].apply(lambda x: x + " Peak area")

    # Melt this into a big table
    pivot_df = input_df.pivot(index='row ID', columns='File', values='Area')

    # Fill NA
    pivot_df = pivot_df.fillna(0)

    # Add in the mz and rt information via a merge to pivot_df
    mz_rt_df = input_df[["row ID", "Mass", "RT"]]
    mz_rt_df = mz_rt_df.drop_duplicates(subset="row ID")
    pivot_df = pd.merge(pivot_df, mz_rt_df, on="row ID")

    # Renaming the columns, m/z
    pivot_df = pivot_df.rename(columns={"Mass": "row m/z", "RT": "row retention time"})

    # Saving the file
    pivot_df.to_csv(output_filename, sep=",", index=False)

    return pivot_df

def _convert_compoundgroups_feature_csv(input_filename, output_filename):
    input_df = pd.read_csv(input_filename, sep='\t')

    new_table_df = pd.DataFrame()

    # We need to turn this table into a proper table
    new_table_df["row ID"] = input_df["Group"]
    new_table_df["row m/z"] = input_df["Mass (avg)"]
    new_table_df["row retention time"] = input_df["RT (avg)"]

    # Figuring out filenames, its all columns after 38
    for column in input_df.columns[38:]:
        new_table_df[column + " Peak area"] = input_df[column]        
    
    return new_table_df

def convert_to_feature_csv(input_filename, output_filename):
    return _convert_compoundgroups_feature_csv(input_filename, output_filename)
    #return _convert_details_feature_csv(input_filename, output_filename)
    

def convert_mgf(input_filenames, compound_feature_table_df, output_mgf):
    compound_map = defaultdict(list)

    for input_filename in input_filenames:
        print(input_filename)

        # reading the mgf with pyteomics
        reader = mgf.read(input_filename)
        for spectrum in reader:
            compound_nr = int(spectrum["params"]["compound_nr"])
            compound_map[compound_nr].append(spectrum)

    output_spectrum_list = []

    # Now we can look at all the features and map to the MS2 spectra
    for feature in compound_feature_table_df.to_dict(orient="records"):
        compound_nr = int(feature["row ID"])

        all_spectra_list = compound_map[compound_nr]
        print(compound_nr, len(all_spectra_list))
        
        if len(all_spectra_list) == 0:
            continue

        min_mz_delta = 10000
        min_spectrum = None
        feature_mass = feature["row m/z"]
        feature_id = feature["row ID"]
        for spectrum in all_spectra_list:
            # Lets choose the one with the closest m/z to the feature id

            # Calculate the delta
            delta = feature_mass - float(spectrum["params"]["pepmass"][0])

            if abs(delta) < min_mz_delta:
                min_mz_delta = abs(delta)
                min_spectrum = spectrum
            
        # Lets add this to the output spectrum list
        # Creating peaks by zipping
        peaks = list(zip(min_spectrum["m/z array"], min_spectrum["intensity array"]))

        new_spectrum = ming_spectrum_library.Spectrum(output_mgf, int(feature_id), int(feature_id), peaks, min_spectrum["params"]["pepmass"][0], 
                                                    min_spectrum["params"]["charge"][0], 2)
        
        output_spectrum_list.append(new_spectrum)

    # Sort list by scan in ascending order
    output_spectrum_list = sorted(output_spectrum_list, key=lambda x: x.scan)

    # Writing the output
    spectrum_collection = ming_spectrum_library.SpectrumCollection("")
    spectrum_collection.spectrum_list = output_spectrum_list
    spectrum_collection.make_scans_contiguous()
    spectrum_collection.save_to_mgf(open(output_mgf, "w"), renumber_scans=False)    

    return spectrum_collection