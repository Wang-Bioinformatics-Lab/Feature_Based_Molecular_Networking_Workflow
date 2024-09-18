#!/usr/bin/python


import sys
import getopt
import os
import json
import argparse
import molecular_network_filtering_library
from collections import defaultdict

def main():
    parser = argparse.ArgumentParser(description='Creating Clustering Info Summary')
    parser.add_argument('networking_pairs_results_file', help='networking_pairs_results_file')
    parser.add_argument('networking_pairs_results_file_filtered', help='networking_pairs_results_file_filtered')
    parser.add_argument('networking_pairs_results_file_filtered_classic_output', help='networking_pairs_results_file_filtered_classic_output')
    parser.add_argument('--top_k_val', help='top_k_val', default=10, type=int)
    parser.add_argument('--max_component_size', help='max_component_size', default=100, type=int)

    args = parser.parse_args()

    top_k_val = args.top_k_val
    max_component_size = args.max_component_size

    G = molecular_network_filtering_library.loading_network(args.networking_pairs_results_file, hasHeaders=True)
    if G == None:
        exit(0)

    molecular_network_filtering_library.filter_top_k(G, top_k_val)
    molecular_network_filtering_library.filter_component(G, max_component_size)
    molecular_network_filtering_library.output_graph_with_headers(G, args.networking_pairs_results_file_filtered)

    #molecular_network_filtering_library.output_graph(G, args.networking_pairs_results_file_filtered_classic_output)


if __name__ == "__main__":
    main()
