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
os.chdir(PROJECT_PATH + "benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'networks', 'load_networks'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbations', 'load_perturbations'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'benchmarking', 'src'))) 
import evaluator
import load_networks
import load_perturbations #not used yet, but soon
importlib.reload(evaluator)
importlib.reload(load_networks)
os.environ["GRN_PATH"]           = PROJECT_PATH + "networks/networks"
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbations/perturbations"

# Parse experiment metadata
with open(os.path.join("experiments", args.experiment_name, "metadata.json")) as f:
    metadata = json.load(f)
assert args.experiment_name == metadata["unique_id"]
outputs = os.path.join("experiments", args.experiment_name, "outputs")
print("\n\n Running experiment " + args.experiment_name + ":\n")
print(yaml.dump(metadata))

# If metadata refers to another experiment, go find code and missing metadata there.
# Otherwise, get the experiment-specific code from here.
if "refers_to" in metadata.keys():
    with open(os.path.join("experiments", metadata["refers_to"], "metadata.json")) as f:
        other_metadata = json.load(f)
    sys.path.append(os.path.expanduser(os.path.join('experiments', other_metadata["unique_id"]))) 
    for key in other_metadata.keys():
        if key not in metadata.keys():
            metadata[key] = other_metadata[key]
else:
    sys.path.append(os.path.expanduser(os.path.join('experiments', args.experiment_name))) 

import this_experiment
importlib.reload(this_experiment)

# Log start time before starting any "real work" involving network structures or pert data
with open(os.path.join(outputs, "..", "start_time.txt"), 'w') as f:
    f.write(str(datetime.datetime.now()) + "\n")

# load networks
def get_subnets(netName):
    print("Getting network '" + netName + "'")
    gc.collect()
    if "random" in netName:
        return {"": evaluator.makeRandomNetwork(density = float(netName[6:])) }
    elif "dense" in netName:
        return {"": evaluator.makeRandomNetwork(density = 1.0)}
    else:
        return {
            subnet: evaluator.networkEdgesToMatrix(load_networks.load_grn_by_subnetwork(netName, subnet))
            for subnet in load_networks.list_subnetworks(netName)
        }


networks = {}
for netName in metadata["network_datasets"]:
    all_subnets = get_subnets(netName)
    for subnet_name in all_subnets.keys():
        new_key = netName + " " + subnet_name if not subnet_name == "" else netName 
        networks[new_key] = all_subnets[subnet_name]
    
# load & split perturbation data
print("Loading & splitting perturbation data.")
perturbed_expression_data = sc.read_h5ad(os.path.join(os.environ["PERTURBATION_PATH"], metadata["perturbation_dataset"], "test.h5ad"))
if args.test_mode:
    perturbed_expression_data = evaluator.downsample(perturbed_expression_data, proportion = 0.2, proportion_genes = 0.01)
allowedRegulators = set.union(*[set(networks[key].columns) for key in networks])
perturbed_expression_data_train, perturbed_expression_data_heldout, perturbationsToPredict = \
    evaluator.splitData(perturbed_expression_data, allowedRegulators, minTestSetSize=5 if args.test_mode else 250)

# Delete any existing models
if args.amount_to_do == "models":
    [os.unlink(os.path.join(outputs, model_file)) for model_file in os.listdir(outputs) if re.search(".celloracle.oracle", model_file)]

# Train or retrieve models (this is memoized)
if args.amount_to_do in {"models", "evaluations"}:
    # Experiment-specific code goes elsewhere
    print("Running experiment")
    experiments, predictions, target_genes, other = this_experiment.run(
        train_data = perturbed_expression_data_train, 
        test_data  = perturbed_expression_data_heldout,
        perturbationsToPredict = perturbationsToPredict,
        networks = networks, 
        outputs = outputs
    )    
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
    #TODO: save predictions?
    stripchartMainFig.figure.savefig(      os.path.join(outputs, "stripchart.pdf"))
    meanSEPlot["spearman"].figure.savefig( os.path.join(outputs, "meanSEPlot.pdf"))
    evaluationResults.to_parquet(          os.path.join(outputs, "networksExperimentEvaluation.parquet"))
    experiments.to_csv(                    os.path.join(outputs, "experiments.csv"))
else: 
    evaluationResults = pd.read_parquet(os.path.join(outputs, "networksExperimentEvaluation.parquet"))

# TODO: remake plots here
if args.amount_to_do in {"plots"}:
    raise NotImplementedError("Remaking just the plots is not yet supported, sorry.")

this_experiment.plot(evaluationResults, outputs)

with open(os.path.join(outputs, "..", "finish_time.txt"), 'w') as f:
    f.write(str(datetime.datetime.now()) + "\n")
