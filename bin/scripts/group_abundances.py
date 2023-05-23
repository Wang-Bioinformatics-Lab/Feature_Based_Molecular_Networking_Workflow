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
    all_attributes = [x for x in metadata_df.columns if x.startswith("ATTRIBUTE_")]

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
    
    # DEBUGGING
    return cluster_summary_df
    
    # Cleaning metadata
    metadata_df["filename"] = metadata_df["filename"].apply(lambda x: os.path.basename(x))

    # Getting filename columns
    filename_columns = [x for x in feature_table_df.columns if x.endswith("Peak area")]

    # Getting the attributes
    all_attributes, all_attribute_groups = create_attribute_group_list(metadata_df)

    for cluster in tqdm(cluster_summary_df.to_dict(orient="records")):
        # Getting the cluster index
        cluster_index = cluster["cluster index"]

        # filter the dataframe
        filtered_cluster_df = feature_table_df[feature_table_df["row ID"] == cluster_index]
        filtered_cluster_df = filtered_cluster_df[filename_columns]

        # transposing the dataframe
        filtered_cluster_df = filtered_cluster_df.transpose()

        # Stripping of Peak area
        filtered_cluster_df["filename"] = [x.replace(" Peak area", "") for x in filtered_cluster_df.index]

        # Setting the area
        filtered_cluster_df["area"] = filtered_cluster_df[0]

        print(metadata_df)

        # Merging in metadata
        filtered_cluster_df = filtered_cluster_df.merge(metadata_df, on="filename", how="inner")

        # Attribute_group
        for attribute_group in all_attribute_groups:
            print(attribute_group)
            attribute = attribute_group["attribute"]
            group_name = attribute_group["group"]

            # Calculating the group mean
            # TODO: Finish This

            cluster[group_column] = group_count

        print(filtered_cluster_df)

        exit(0)

        break


    # Getting all the filenames
    #print(metadata_df)

    print("XXX")


    # Cleaning the filenames
    clusterinfo_df["#Filename"] = clusterinfo_df["#Filename"].apply(lambda x: os.path.basename(x))
    metadata_df["filename"] = metadata_df["filename"].apply(lambda x: os.path.basename(x))

    # Folding in the metadata into clusterinfo
    clusterinfo_df = clusterinfo_df.merge(metadata_df, left_on="#Filename", right_on="filename", how="left")

    # First lets group the cluster info by cluster index
    grouped_clusterinfo_df = clusterinfo_df.groupby("#ClusterIdx")

    # loop through all the clusters
    cluster_summary_list = clustersummary_df.to_dict(orient="records")

    

    for cluster in tqdm(cluster_summary_list):
        # filter for the grouped list
        cluster_index = cluster["cluster index"]
        clusterinfo_per_group_df = grouped_clusterinfo_df.get_group(cluster_index)

        # TODO: We can likely speed this up with pandas operations

        # Attribute_group
        for attribute_group in all_attribute_groups:
            #print(attribute_group)

            # filtering the data
            group_count = len(clusterinfo_per_group_df[clusterinfo_per_group_df[attribute_group["attribute"]] == attribute_group["group"]])
            group_column = "{}:GNPSGROUP:{}".format(attribute_group["attribute"], attribute_group["group"])

            cluster[group_column] = group_count
    
        # Adding the cluster information for which group membership it is a part of
        for attribute in all_attributes:
            # Finding all groups in the attribute
            all_groups = set(clusterinfo_per_group_df[attribute])
            cluster[attribute] = ",".join(all_groups)

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
