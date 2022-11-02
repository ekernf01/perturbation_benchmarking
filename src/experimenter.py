
import os
import re
import matplotlib.colors as colors
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
import pandas as pd
import scanpy as sc
import gc
import json
import yaml
import argparse
from argparse import Namespace
import datetime

# User input: name of experiment and whether to fully rerun or just remake plots. 
parser = argparse.ArgumentParser("experimenter")
parser.add_argument("--experiment_name", help="Unique id for the experiment.", type=str)
parser.add_argument("--test_mode",       help="If true, use a small subset of the perturbation data.",     default = False, action = "store_true")
parser.add_argument(
    "--amount_to_do",
    choices = ["plots", "evaluations", "models"],
    help="""
    The code makes models, evaluations, and plots, in that order. It saves the models and the evaluations. 
    To do just plots, using saved evaluations and models, specify "plots".
    To do plots and evaluations using saved models, specify "evaluations".
    To do everything, specify "models". 
    """
)
args = parser.parse_args()
print("args to experimenter.py:")
print(args)

# For interactive sessions
if args.experiment_name is None:
    args = Namespace(**{
        "experiment_name":'test',
        "test_mode":True,
        "amount_to_do": "models"
    })

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
outputs = os.path.join("experiments", args.experiment_name, "outputs")

# Parse experiment metadata
def validate_metadata():
    with open(os.path.join("experiments", args.experiment_name, "metadata.json")) as f:
        metadata = json.load(f)
    if not metadata["is_active"]:
        raise ValueError("This experiment is marked as inactive. If you really want to run it, edit its metadata.json.")
    print("\n\n Running experiment " + args.experiment_name + ":\n")
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
        code_location = os.path.expanduser(os.path.join('experiments', args.experiment_name))

    # network handling is complex; add some default behavior to reduce metadata boilerplate
    for netName in metadata["network_datasets"].keys():
        if not "subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["subnets"] = ["all"]
        if not "do_aggregate_subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["do_aggregate_subnets"] = False
    
    # Check basics on metadata
    assert args.experiment_name == metadata["unique_id"], "Experiment is labeled right"
    assert metadata["perturbation_dataset"] in set(load_perturbations.load_perturbation_metadata().query("is_ready=='yes'")["name"]), "perturbation data exist as named"
    for netName in metadata["network_datasets"].keys():
        assert netName in set(load_networks.load_grn_metadata()["name"]).union({"dense"}) or "random" in netName, "Networks exist as named"
        assert "subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"
        assert "do_aggregate_subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"

    print("\nFully parsed metadata:\n")
    print(yaml.dump(metadata))

    return metadata, code_location

# Load experiment code
metadata, code_location = validate_metadata()
sys.path.append(code_location)
import this_experiment
importlib.reload(this_experiment)

# Log start time before starting any "real work" involving network structures or pert data
with open(os.path.join(outputs, "..", "start_time.txt"), 'w') as f:
    f.write(str(datetime.datetime.now()) + "\n")

# Delete any existing models if models are to be re-fitted
if args.amount_to_do == "models":
    [os.unlink(os.path.join(outputs, model_file)) for model_file in os.listdir(outputs) if re.search(".celloracle.oracle", model_file)]

# load networks
# TODO: move this to the network collection loader module
def get_subnets(netName:str, subnets:list, test_mode: bool = args.test_mode) -> dict:
    """Get gene regulatory networks.

    Args:
        netName (str): Name of network to pull from collection, or "dense" or e.g. "random0.123" for random with density 12.3%. 
        subnets (list, optional): List of cell type- or tissue-specific subnetworks to include. 
        test_mode (bool, optional): Lighten the load during testing. Defaults to args.test_mode.

    Returns:
        dict: A dict containing base GRN's in the format expected by CO.
    """
    print("Getting network '" + netName + "'")
    gc.collect()
    if "random" in netName:
        return { "": evaluator.pivotNetworkWideToLong( evaluator.makeRandomNetwork( density = float( netName[6:] ) ) ) }
    elif "dense" in netName:
        return { "": evaluator.pivotNetworkWideToLong( evaluator.makeRandomNetwork( density = 1.0 ) ) }
    else:
        if subnets[0]=="all":
            subnets = load_networks.list_subnetworks(netName)
        if test_mode:
            subnets = subnets[0:1]
        return {
            subnet: load_networks.load_grn_by_subnetwork(netName, subnet)
            for subnet in subnets
        }

if args.amount_to_do in {"models", "evaluations"}:
    # Get networks
    networks = {}
    for netName in metadata["network_datasets"].keys():
        all_subnets = get_subnets(netName, subnets = metadata["network_datasets"][netName]["subnets"] )
        if metadata["network_datasets"][netName]["do_aggregate_subnets"]:
            raise NotImplementedError("Sorry, still haven't implemented unions of network edge sets. For now, set do_aggregate_subnets to False.")
        else:
            for subnet_name in all_subnets.keys():
                new_key = netName + " " + subnet_name if not subnet_name == "" else netName 
                networks[new_key] = all_subnets[subnet_name]

    # load & split perturbation data
    print("Loading & splitting perturbation data.")
    perturbed_expression_data = load_perturbations.load_perturbation(metadata["perturbation_dataset"])
    if args.test_mode:
        perturbed_expression_data = evaluator.downsample(perturbed_expression_data, proportion = 0.2, proportion_genes = 0.01)
    if any(networks):
        allowedRegulators = set.union(*[set(networks[key]["regulator"]) for key in networks])
    else:
        allowedRegulators = perturbed_expression_data.var_names
    perturbed_expression_data_train, perturbed_expression_data_heldout = \
        evaluator.splitData(perturbed_expression_data, allowedRegulators, minTestSetSize=5 if args.test_mode else 250)

    # Experiment-specific code goes elsewhere
    print("Running experiment")
    pp = [
            (r[1][0], r[1][1]) 
            for r in perturbed_expression_data_heldout.obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
        ]
    experiments, predictions, other = this_experiment.run(
        train_data = perturbed_expression_data_train, 
        test_data  = perturbed_expression_data_heldout,
        perturbationsToPredict = pp,
        networks = networks, 
        outputs = outputs
    )      
    experiments.to_csv( os.path.join(outputs, "experiments.csv") )
    print("Evaluating predictions.")
    evaluationResults, stripchartMainFig, meanSEPlot = evaluator.evaluateCausalModel(
        heldout = perturbed_expression_data_heldout, 
        predictions = predictions, 
        baseline = perturbed_expression_data_train[perturbed_expression_data_train.obs["is_control"], :],
        experiments = experiments, 
        outputs = outputs, 
        factor_varied = metadata["factor_varied"], 
        default_level = metadata["default_level"], 
        classifier = None
    )
    evaluationResults.to_parquet(          os.path.join(outputs, "networksExperimentEvaluation.parquet"))

if args.amount_to_do in {"plots", "models", "evaluations"}:
    evaluationResults = pd.read_parquet(os.path.join(outputs, "networksExperimentEvaluation.parquet"))
    evaluator.makeMainPlots(
        evaluationResults, 
        factor_varied = metadata["factor_varied"],
        outputs = outputs, 
        default_level=metadata["default_level"],
    )


this_experiment.plot(evaluationResults, outputs)

with open(os.path.join(outputs, "..", "finish_time.txt"), 'w') as f:
    f.write(str(datetime.datetime.now()) + "\n")
