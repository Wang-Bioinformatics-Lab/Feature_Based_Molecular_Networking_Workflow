#!/usr/bin/python


import sys
import getopt
import os
import argparse
import statistics
import pandas as pd
from tqdm import tqdm
from collections import defaultdict


def create_attribute_group_list(metadata_df):
    # Determining all the groups we want to calculate over columsn with prefix ATTRIBUTE_
    all_attributes = [x for x in metadata_df.columns if x.upper().startswith("ATTRIBUTE_")]

    all_attribute_groups = []

    for attribute in all_attributes:
        # Getting the groups in each attribute and creating a list of them
        attribute_groups = metadata_df[attribute].unique().tolist()

        # Creating a dictionary for each attribute group
        for attribute_group in attribute_groups:
            attribute_group_dict = {}
            attribute_group_dict["attribute"] = attribute
            attribute_group_dict["group"] = attribute_group
            all_attribute_groups.append(attribute_group_dict)

    return all_attributes, all_attribute_groups


# This function calculates all the group counds for all the relevant columns
def calculate_groups_metadata(feature_table_df, metadata_df):
    # Creating clustersummary
    cluster_summary_df = pd.DataFrame()
    cluster_summary_df["cluster index"] = feature_table_df["row ID"]
    cluster_summary_df["parent mass"] = feature_table_df["row m/z"]

    if len(metadata_df) == 0:
        return cluster_summary_df
    
    if "filename" not in metadata_df.columns:
        raise Exception("Metadata does not contain filename column")
    
    # Cleaning metadata
    metadata_df["filename"] = metadata_df["filename"].apply(lambda x: os.path.basename(x))

    # Getting filename columns
    filename_columns = [x for x in feature_table_df.columns if x.endswith("Peak area")]

    # Getting the attributes
    all_attributes, all_attribute_groups = create_attribute_group_list(metadata_df)

    # Making the quant table tall
    tall_feature_table_df = pd.melt(feature_table_df, id_vars=["row ID", "row m/z", "row retention time"], value_vars=filename_columns, var_name="filename", value_name="area")
    tall_feature_table_df["filename"] = tall_feature_table_df["filename"].apply(lambda x: x.replace("Peak area", "").rstrip())

    # Merging in metadata
    tall_feature_table_df = tall_feature_table_df.merge(metadata_df, on="filename", how="inner")


    # Doing the actual calculations
    cluster_summary_list = cluster_summary_df.to_dict(orient="records")

    for cluster in tqdm(cluster_summary_list):
        # Getting the cluster index
        cluster_index = cluster["cluster index"]

        # filter the dataframe
        filtered_cluster_df = tall_feature_table_df[tall_feature_table_df["row ID"] == cluster_index]

        # Attribute_group
        for attribute_group in all_attribute_groups:
            attribute = attribute_group["attribute"]
            group_name = attribute_group["group"]

            try:
                # filtering the data
                group_data_df = filtered_cluster_df[filtered_cluster_df[attribute] == group_name]

                # Calcualting the average area
                area_average = group_data_df["area"].mean()
            except:
                area_average = 0
            
            group_column = "{}:GNPSGROUP:{}".format(attribute, group_name)

            cluster[group_column] = area_average

    return pd.DataFrame(cluster_summary_list)




def main():
    parser = argparse.ArgumentParser(description='Creates enriched cluster info summary')
    parser.add_argument('input_featuretable', help='input_featuretable')
    parser.add_argument('input_metadata', help='input_group_mapping_filename')
    parser.add_argument('output_clusterinfosummary_filename', help='output_clusterinfosummary_filename')
    args = parser.parse_args()

    # Loading Data
    feature_table_df = pd.read_csv(args.input_featuretable, sep=",")

    try:
        metadata_df = pd.read_csv(args.input_metadata, sep="\t")
    except:
        metadata_df = pd.DataFrame()

    print(metadata_df)

    # Enriching metadata group counts
    # TODO :finish this efficiently
    clustersummary_df = calculate_groups_metadata(feature_table_df, metadata_df)

    # Writing out the file
    clustersummary_df.to_csv(args.output_clusterinfosummary_filename, sep="\t", index=False)

    exit(0)




if __name__ == "__main__":
    main()
