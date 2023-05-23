#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mingxun wang
@purpose: to convert the Agilent file into a diserable format
"""
import pandas as pd
import sys

def convert_to_feature_csv(input_filename, output_filename):

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