
import os
import gc
import json
import yaml
import gc
import pandas as pd
import anndata
from itertools import product
# Deal with various file paths specific to this project
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
import evaluator
import ggrn
import load_networks
import load_perturbations
importlib.reload(evaluator)
importlib.reload(load_networks)
importlib.reload(load_perturbations)
os.environ["GRN_PATH"]           = PROJECT_PATH + "network_collection/networks"
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbation_data/perturbations"

# Parse experiment metadata
def validate_metadata(
    experiment_name, 
    permissive = False
):
    with open(os.path.join("experiments", experiment_name, "metadata.json")) as f:
        metadata = json.load(f)
    if not permissive and not metadata["is_active"]:
        raise ValueError("This experiment is marked as inactive. If you really want to run it, edit its metadata.json.")
    print("\n\nRaw metadata for experiment " + experiment_name + ":\n")
    print(yaml.dump(metadata))

    # If metadata refers to another experiment, go find missing metadata there.
    if "refers_to" in metadata.keys():
        with open(os.path.join("experiments", metadata["refers_to"], "metadata.json")) as f:
            other_metadata = json.load(f)
            assert other_metadata["is_active"], "Referring to an inactive experiment is not allowed."
        for key in other_metadata.keys():
            if key not in metadata.keys():
                metadata[key] = other_metadata[key]
    else:
        metadata["refers_to"] = None

    # Set defaults (None often defers to downstream code)
    defaults = {
        "pruning_parameter": None, 
        "pruning_strategy": None,
        "network_prior": None,
        "desired_heldout_fraction": 0.5,
        "type_of_split": "interventional",
        "regression_method": "RidgeCV",
        "time_strategy": "steady_state",
        "kwargs": None,
    }
    for k in defaults:
        if not k in metadata:
            metadata[k] = defaults[k]

    # network handling is complex; add some default behavior to reduce metadata boilerplate
    for netName in metadata["network_datasets"].keys():
        if not "subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["subnets"] = ["all"]
        if not "do_aggregate_subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["do_aggregate_subnets"] = False
    
    # Check all keys
    missing = [k for k in (
        # Experiment info
        "unique_id",
        "nickname",
        "readme",
        "question",
        "is_active",
        "factor_varied",    
        "color_by",
        "facet_by",
        # Data and preprocessing
        "network_datasets",
        "perturbation_dataset",
        "merge_replicates",
        "desired_heldout_fraction",
        "type_of_split",
        # Modeling decisions
        "pruning_parameter", 
        "pruning_strategy",
        "network_prior",
        "regression_method",
        "time_strategy",
    ) if k not in metadata.keys()]
    assert len(missing)==0, f"Metadata is missing some required keys: {' '.join(missing)}"
    
    # Check a few of the values
    assert experiment_name == metadata["unique_id"], "Experiment is labeled right"
    if not permissive:
        assert metadata["perturbation_dataset"] in set(load_perturbations.load_perturbation_metadata().query("is_ready=='yes'")["name"]), "perturbation data exist as named"
        for netName in metadata["network_datasets"].keys():
            assert netName in set(load_networks.load_grn_metadata()["name"]).union({"dense", "empty"}) or "random" in netName, "Networks exist as named"
            assert "subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"
            assert "do_aggregate_subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"

    print("\nFully parsed metadata:\n")
    print(yaml.dump(metadata))

    return metadata

def lay_out_runs(
  networks: dict, 
  metadata: dict,
) -> pd.DataFrame:
    """Lay out the specific runs done in this experiment.

    Args:
    networks (dict): dict with string keys and LightNetwork values
    outputs (str): folder name to save results in
    metadata (dict): metadata for this Experiment, from metadata.json. See this repo's global README.

    Returns:
        pd.DataFrame: metadata on the different conditions in this experiment

    """
    metadata = metadata.copy() # We're gonna mangle it. :)
    # This will fuck up the cartesian product below.
    del metadata["kwargs"]
    # See experimenter.get_networks() to see how the metadata.json evolves into this thing
    metadata["network"] = list(networks.keys())
    # Code downstream (product) splits strings if you don't do this.
    for k in metadata.keys():
        if type(metadata[k]) != list:
            metadata[k] = [metadata[k]]
    # Combos 
    experiments =  pd.DataFrame(
        [row for row in product(*metadata.values())], 
        columns=metadata.keys()
    )
    # The dense network is represented as an empty network to save space.
    # Recommended usage is to set network_prior="ignore", otherwise the empty network will be taken literally. 
    for i in experiments.index:
        experiments.loc[i, "network_prior"] = \
        "ignore" if experiments.loc[i, "network"] == "dense" else experiments.loc[i, "network_prior"]

    return experiments
  
def do_one_run(
  experiments: pd.DataFrame, 
  i: int, 
  train_data: anndata.AnnData, 
  test_data: anndata.AnnData, 
  networks: dict, 
  outputs: str,
  metadata: dict, 
  ) -> anndata.AnnData:
  """Do one run (fit a GRN model and make predictions) as part of this experiment.

  Args:
      experiments (pd.DataFrame): Output of lay_out_runs
      i (int): A value from the experiments.index
      Other args: see help(lay_out_runs)

  Returns:
      anndata.AnnData: Predicted expression
  """
  grn = ggrn.GRN(
    train=train_data, 
    network=networks[experiments.loc[i,'network']]
  )
  grn.extract_tf_activity(method = "tf_rna")
  grn.fit(
      method                               = experiments.loc[i,"regression_method"], 
      cell_type_labels                     = None,
      cell_type_sharing_strategy           = "identical",
      network_prior                        = experiments.loc[i,"network_prior"],
      pruning_strategy                     = experiments.loc[i,"pruning_strategy"],
      pruning_parameter                    = experiments.loc[i,"pruning_parameter"],
      projection                           = "none",  
      time_strategy                        = experiments.loc[i,"time_strategy"],
      kwargs                               = metadata["kwargs"],
  )
  return grn

# TODO: move this to the network collection loader module?
def get_subnets(netName:str, subnets:list, test_mode, target_genes = None, do_aggregate_subnets = False) -> dict:
    """Get gene regulatory networks for an experiment.

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
    elif "empty" == netName or "dense" == netName:
        networks = { 
            netName: load_networks.LightNetwork(df=pd.DataFrame(index=[], columns=["regulator", "target", "weight"]))
        }
        if "dense"==netName:
            print("WARNING: for 'dense' network, returning an empty network. In GRN.fit(), use network_prior='ignore'. ")
    else:            
        networks = {}
        if do_aggregate_subnets:
            new_key = netName 
            if subnets[0]=="all":
                networks[new_key] = load_networks.LightNetwork(netName)
            else:
                networks[new_key] = load_networks.LightNetwork(netName, subnets)
        else:
            for subnet_name in subnets:
                new_key = netName + " " + subnet_name
                if subnets[0]=="all":
                    networks[new_key] = load_networks.LightNetwork(netName)
                else:
                    networks[new_key] = load_networks.LightNetwork(netName, [subnet_name])
    return networks


