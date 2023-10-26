#!/usr/bin/python


import sys
import getopt
import os
import molecular_network_filtering_library
import networkx as nx
import argparse


def convert_network(G):
    # This helps to reformat from the standard GNPS library to be able to display in GNPS2
    new_G = nx.Graph()

    all_nodes = list(G.nodes())

    group_columns = G.nodes[all_nodes[0]].keys()
    group_columns = [key for key in group_columns if "GNPSGROUP" in key]

    node_to_component = {}

    for node in all_nodes:
        new_G.add_node(node)

        # Fixing Group Names
        for column in group_columns:
            if column.upper().startswith("ATTRIBUTE"):
                new_key = column.upper()
            else:
                new_key = "ATTRIBUTE_GNPS:{}".format(column.replace("GNPSGROUP:", ""))
            new_G.nodes[node][new_key] = float(G.nodes[node][column])

        # Fixing node attributes
        new_G.nodes[node]["mz"] = float("{:.4f}".format(float(G.nodes[node]["parent mass"])))
        new_G.nodes[node]["rt"] = float("{:.2f}".format(float(G.nodes[node]["RTMean"])))
        new_G.nodes[node]["rt_min"] = float("{:.2f}".format(float(G.nodes[node]["RTMean"])))
        new_G.nodes[node]["charge"] = G.nodes[node]["charge"]
        new_G.nodes[node]["component"] = G.nodes[node]["component"]

        node_to_component[node] = G.nodes[node]["component"]
        
        
        if "Compound_Name" in G.nodes[node]:
            new_G.nodes[node]["library_compound_name"] = G.nodes[node]["Compound_Name"]
            new_G.nodes[node]["library_SMILES"] = G.nodes[node]["Smiles"]
            new_G.nodes[node]["library_InChI"] = G.nodes[node]["INCHI"]

    # Fixing Edges
    for node1, node2, data in G.edges.data():

        if node1 != node2:
            new_G.add_edge(node1, node2)

            # Get edge attributes
            new_G[node1][node2]["deltamz"] = float(data["mass_difference"])
            new_G[node1][node2]["deltamz_int"] = int(float(data["mass_difference"]))
            new_G[node1][node2]["score"] = float(data["cosine_score"])
            new_G[node1][node2]["matched_peaks"] = "0"
            new_G[node1][node2]["scan1"] = node1
            new_G[node1][node2]["scan2"] = node2
            new_G[node1][node2]["component"] = node_to_component[node2]

    return new_G


def main():
    parser = argparse.ArgumentParser(description='Creating Clustering Info Summary')
    parser.add_argument('input_clusterinfo_summary', help='input_clusterinfo_summary')
    parser.add_argument('input_pairs', help='input_pairs')
    parser.add_argument('input_library_matches', help='input_library_matches')
    parser.add_argument('output_graphml', help='output_graphml')

    args = parser.parse_args()

    #Doing other filtering
    G = molecular_network_filtering_library.loading_network(args.input_pairs, hasHeaders=True)
    molecular_network_filtering_library.add_clusterinfo_summary_to_graph(G, args.input_clusterinfo_summary)
    molecular_network_filtering_library.add_library_search_results_to_graph(G, args.input_library_matches)

    # Cleaning up network when the clusterinfo summary is not present 

    # Reformatting
    G = convert_network(G)

    nx.write_graphml(G, args.output_graphml, infer_numeric_types=True)




if __name__ == "__main__":
    main()
