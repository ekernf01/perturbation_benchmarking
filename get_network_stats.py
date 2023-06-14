# This script computes the in-degree of each target and out-degree of each regulator across all our pre-compiled GRN's.
# For collections with multiple tissues, degree is the average over tissues, so it may be a fraction.

import os
import pandas as pd
import numpy as np
import scanpy as sc
# Deal with various file paths specific to this project
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import load_networks
importlib.reload(load_networks)
os.environ["GRN_PATH"]           = PROJECT_PATH + "network_collection/networks"
degree_info = {}
for network_name in load_networks.load_grn_metadata()["name"]:
    print(network_name)
    num_subnetworks = len(load_networks.list_subnetworks(network_name))
    try:
        network_edges = load_networks.load_grn_all_subnetworks(network_name)
    except Exception as e:
        print(f"Loading network {network_name} failed will error: \n{repr(e)}")
        continue
    outdegrees = network_edges["regulator"].value_counts()    
    indegrees = network_edges["target"].value_counts()
    degree_info[network_name] = pd.merge(outdegrees, indegrees, left_index=True, right_index=True, how = "outer").fillna(0)
    degree_info[network_name].columns = ["out-degree", "in-degree"]
    degree_info[network_name] = degree_info[network_name]/num_subnetworks 
    degree_info[network_name]["network"] = network_name
degree_info = pd.concat(degree_info.values())
degree_info.to_csv("../accessory_data/degree_info.csv")
