
import os
import gc
import json
import yaml
import gc
# Deal with various file paths specific to this project
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
import evaluator
import load_networks
import load_perturbations
importlib.reload(evaluator)
importlib.reload(load_networks)
importlib.reload(load_perturbations)
os.environ["GRN_PATH"]           = PROJECT_PATH + "network_collection/networks"
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbation_data/perturbations"

# Parse experiment metadata
def validate_metadata(experiment_name, permissive = False):
    with open(os.path.join("experiments", experiment_name, "metadata.json")) as f:
        metadata = json.load(f)
    if not permissive and not metadata["is_active"]:
        raise ValueError("This experiment is marked as inactive. If you really want to run it, edit its metadata.json.")
    print("\n\n Raw metadata experiment " + experiment_name + ":\n")
    print(yaml.dump(metadata))

    # If metadata refers to another experiment, go find code and missing metadata there.
    # Otherwise, get the experiment-specific code from here.
    if "refers_to" in metadata.keys():
        with open(os.path.join("experiments", metadata["refers_to"], "metadata.json")) as f:
            other_metadata = json.load(f)
        code_location = os.path.expanduser(os.path.join('experiments', other_metadata["unique_id"]))
        for key in other_metadata.keys():
            if key not in metadata.keys():
                metadata[key] = other_metadata[key]
    else:
        code_location = os.path.expanduser(os.path.join('experiments', experiment_name))
        metadata["refers_to"] = None

    # network handling is complex; add some default behavior to reduce metadata boilerplate
    for netName in metadata["network_datasets"].keys():
        if not "subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["subnets"] = ["all"]
        if not "do_aggregate_subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["do_aggregate_subnets"] = False
    
    # Check basics on metadata
    assert experiment_name == metadata["unique_id"], "Experiment is labeled right"
    if not permissive:
        assert metadata["perturbation_dataset"] in set(load_perturbations.load_perturbation_metadata().query("is_ready=='yes'")["name"]), "perturbation data exist as named"
        for netName in metadata["network_datasets"].keys():
            assert netName in set(load_networks.load_grn_metadata()["name"]).union({"dense"}) or "random" in netName, "Networks exist as named"
            assert "subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"
            assert "do_aggregate_subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"

    print("\nFully parsed metadata:\n")
    print(yaml.dump(metadata))

    return metadata, code_location


# TODO: move this to the network collection loader module?
def get_subnets(netName:str, subnets:list, test_mode, target_genes = None, do_aggregate_subnets = False) -> dict:
    """Get gene regulatory networks.

    Args:
        netName (str): Name of network to pull from collection, or "dense" or e.g. "random0.123" for random with density 12.3%. 
        subnets (list, optional): List of cell type- or tissue-specific subnetworks to include. 
        test_mode (bool, optional): Lighten the load during testing. Defaults to args.test_mode.
        do_aggregate_subnets (bool, optional): If True, return has just one network named netName. If False,
            then returned dict has many separate networks named like netName + " " + subnet_name.

    Returns:
        dict: A dict containing base GRN's as LightNetwork objects (see the docs in the load_networks module in the networks collection.)
    """
    print("Getting network '" + netName + "'")
    gc.collect()
    if "random" in netName:
        networks = { 
            netName: load_networks.LightNetwork(
                df = evaluator.pivotNetworkWideToLong( 
                    evaluator.makeRandomNetwork( target_genes = target_genes, density = float( netName[6:] ) ) 
                ) 
            )
        }
    else:
        if do_aggregate_subnets:
            if subnets[0]=="all":
                networks = load_networks.LightNetwork(netName)
            else:
                networks = load_networks.LightNetwork(netName, subnets)
        else:
            networks = {}
            for subnet_name in subnets:
                new_key = netName + " " + subnet_name
                networks[new_key] = load_networks.LightNetwork(netName, subnet_name)
    return networks