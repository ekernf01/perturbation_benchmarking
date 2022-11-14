
import os
import numpy as np
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
import gc
# User input: name of experiment and whether to fully rerun or just remake plots. 
parser = argparse.ArgumentParser("experimenter")
parser.add_argument("--experiment_name", help="Unique id for the experiment.", type=str)
parser.add_argument("--test_mode",       help="If true, use a small subset of the perturbation data.",     default = False, action = "store_true")
parser.add_argument(
    "--amount_to_do",
    choices = ["plots", "evaluations", "models", "missing_models"],
    help="""
    The code makes models, evaluations, and plots, in that order. It saves the models and the evaluations. 
    To do just plots, using saved evaluations and models, specify "plots".
    To do plots and evaluations using saved models, specify "evaluations".
    To do everything, specify "models". 
    If it crashes, specify "missing_models" to keep previous progress. 
    """
)
args = parser.parse_args()
print("args to experimenter.py:", flush = True)
print(args)

# For interactive sessions
if args.experiment_name is None:
    args = Namespace(**{
        "experiment_name":'1.4_1',
        "test_mode":False,
        "amount_to_do": "missing_models"
    })

# Deal with various file paths specific to this project
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
import evaluator
import experimenter
import load_perturbations
importlib.reload(evaluator)
importlib.reload(experimenter)
importlib.reload(load_perturbations)
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbation_data/perturbations"
outputs = os.path.join("experiments", args.experiment_name, "outputs")


# Load experiment code
metadata, code_location = experimenter.validate_metadata(experiment_name=args.experiment_name)
sys.path.append(code_location)
import this_experiment
importlib.reload(this_experiment)

print("Starting at " + str(datetime.datetime.now()), flush = True)

if args.amount_to_do in {"models", "missing_models", "evaluations"}:
    perturbed_expression_data = load_perturbations.load_perturbation(metadata["perturbation_dataset"])
    # perturbed_expression_data = evaluator.averageWithinPerturbation(perturbed_expression_data)

    # Get networks
    networks = {}
    for netName in metadata["network_datasets"].keys():
        networks = networks | experimenter.get_subnets(
            netName, 
            subnets = metadata["network_datasets"][netName]["subnets"], 
            target_genes = perturbed_expression_data.var_names, 
            test_mode = args.test_mode, 
            do_aggregate_subnets = metadata["network_datasets"][netName]["do_aggregate_subnets"]
        )

    # load & split perturbation data
    print("Loading & splitting perturbation data.", flush = True)
    try:
        perturbed_expression_data = perturbed_expression_data.to_memory()
    except ValueError: #Object is already in memory.
        pass

    if args.test_mode:
        perturbed_expression_data = evaluator.downsample(perturbed_expression_data, proportion = 0.2, proportion_genes = 0.01)
    if any(networks):
        allowedRegulators = set.union(*[networks[key].get_all_regulators() for key in networks])
    else:
        allowedRegulators = perturbed_expression_data.var_names
    perturbed_expression_data_train, perturbed_expression_data_heldout = \
        evaluator.splitData(perturbed_expression_data, allowedRegulators, desired_heldout_fraction=0.5)
    del perturbed_expression_data
    gc.collect()
    # Experiment-specific code goes elsewhere
    pp = [
            (r[1][0], r[1][1]) 
            for r in perturbed_expression_data_heldout.obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
        ]

if args.amount_to_do in {"models", "missing_models"}:
    print("Running experiment", flush = True)
    experiments = this_experiment.lay_out_runs(
            train_data = perturbed_expression_data_train, 
            test_data  = perturbed_expression_data_heldout,
            perturbationsToPredict = pp,
            networks = networks, 
            outputs = outputs,
            )
    experiments.to_csv( os.path.join(outputs, "experiments.csv") )
    os.makedirs(os.path.join( outputs, "predictions" ), exist_ok=True) 
    for i in experiments.index:
        h5ad = os.path.join( outputs, "predictions", str(i) + ".h5ad" )
        if (args.amount_to_do in {"models"}) or not os.path.isfile(h5ad):
            try:
                os.unlink(h5ad)
            except FileNotFoundError:
                pass
            predictions = this_experiment.do_one_run(
                experiments, 
                i,
                train_data = perturbed_expression_data_train, 
                test_data  = perturbed_expression_data_heldout,
                perturbationsToPredict = pp,
                networks = networks, 
                outputs = outputs,
            )
            print("Saving predictions...", flush = True)
            predictions.write_h5ad( h5ad )
            print("... done.", flush = True)

if args.amount_to_do in {"models", "missing_models", "evaluations"}:
    print("Retrieving saved predictions", flush = True)
    experiments = pd.read_csv( os.path.join(outputs, "experiments.csv") )
    predictions = {i:sc.read_h5ad( os.path.join(outputs, "predictions", str(i) + ".h5ad" ) ) for i in experiments.index}

    evaluationResults, vlnplot = evaluator.evaluateCausalModel(
        heldout = perturbed_expression_data_heldout, 
        predictions = predictions, 
        baseline = perturbed_expression_data_train[perturbed_expression_data_train.obs["is_control"], :],
        experiments = experiments, 
        outputs = outputs, 
        factor_varied = metadata["factor_varied"], 
        default_level = None,
        classifier = None,
    )
    evaluationResults.to_parquet(          os.path.join(outputs, "networksExperimentEvaluation.parquet"))

if args.amount_to_do in {"plots", "models", "missing_models", "evaluations"}:
    evaluationResults = pd.read_parquet(os.path.join(outputs, "networksExperimentEvaluation.parquet"))
    evaluator.makeMainPlots(
        evaluationResults, 
        outputs = outputs, 
        factor_varied = metadata["factor_varied"],
        color_by = metadata["color_by"],
        facet_by=metadata["facet_by"],
    )


this_experiment.plot(evaluationResults, outputs)

print("Experiment done at " + str(datetime.datetime.now()), flush = True)