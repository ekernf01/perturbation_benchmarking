import os
import numpy as np
import pandas as pd
import scanpy as sc
import gc
import argparse
from argparse import Namespace
import datetime

# Access our code
import perturbation_benchmarking_package.evaluator as evaluator
import perturbation_benchmarking_package.experimenter as experimenter
import load_perturbations
import load_networks
import ggrn.api as ggrn

# Access our data collections
load_networks.set_grn_location(
    '../network_collection/networks'
)
load_perturbations.set_data_path(
    '../perturbation_data/perturbations'
)
DEFAULT_HUMAN_TFs = pd.read_csv("../accessory_data/humanTFs.csv")
DEFAULT_HUMAN_TFs = DEFAULT_HUMAN_TFs.loc[DEFAULT_HUMAN_TFs["Is TF?"]=="Yes", "HGNC symbol"]

# User input: name of experiment and whether to fully rerun or just remake plots. 
parser = argparse.ArgumentParser("experimenter")
parser.add_argument("--experiment_name", help="Unique id for the experiment.", type=str)
parser.add_argument("--save_models",     help="If true, save model objects.", default = False, action = "store_true")
parser.add_argument("--save_trainset_predictions", help="If true, make & save predictions of training data.", default = False, action = "store_true")
parser.add_argument('--no_skip_bad_runs', dest='skip_bad_runs', action='store_false', help="Unless this flag is used, keep running when some runs hit errors.")
parser.set_defaults(feature=True)
parser.add_argument(
    "--amount_to_do",
    choices = ["plots", "evaluations", "models", "missing_models"],
    help="""
    The code makes models, evaluations, and plots, in that order. It saves the models and the evaluations. 
    To do just plots, using saved evaluations and models, specify "plots".
    To do plots and evaluations using saved models, specify "evaluations".
    To do everything, specify "models". 
    If it crashes, specify "missing_models" to keep previous progress. 
    To skip certain models (e.g. skip ExtraTrees if low on RAM), manually place 
    an empty results file like 'touch outputs/results/predictions/3.h5ad'.
    """
)
args = parser.parse_args()
print("args to experimenter.py:", flush = True)
print(args)
# Default args to this script for interactive use
if args.experiment_name is None:
    args = Namespace(**{
        "experiment_name": "1.0_0",
        "amount_to_do": "missing_models",
        "save_trainset_predictions": True,
        "save_models": False,
        "skip_bad_runs": False,
    })
# Additional bookkeeping
print("Running experiment", flush = True)
outputs = os.path.join("experiments", args.experiment_name, "outputs")
os.makedirs(outputs, exist_ok=True)
metadata = experimenter.validate_metadata(experiment_name=args.experiment_name)
print("Starting at " + str(datetime.datetime.now()), flush = True)

# Set up the perturbation and network data
perturbed_expression_data, networks, conditions = experimenter.set_up_data_networks_conditions(
    metadata,
    amount_to_do = args.amount_to_do, 
    outputs = outputs,
)
def get_current_data_split(i, verbose = False):
    return experimenter.splitDataWrapper(
        experimenter.filter_genes(perturbed_expression_data, num_genes = conditions.loc[i, "num_genes"], outputs = outputs),
        networks = networks, 
        desired_heldout_fraction = conditions.loc[i, "desired_heldout_fraction"],  
        type_of_split            = conditions.loc[i, "type_of_split"],
        data_split_seed          = conditions.loc[i, "data_split_seed"],
        verbose = verbose,
    )

# Begin conditions
os.makedirs(os.path.join( outputs, "predictions"   ), exist_ok=True) 
os.makedirs(os.path.join( outputs, "fitted_values" ), exist_ok=True) 
for i in conditions.index:
    models      = os.path.join( outputs, "models",        str(i) )
    h5ad        = os.path.join( outputs, "predictions",   str(i) + ".h5ad" )
    h5ad_fitted = os.path.join( outputs, "fitted_values", str(i) + ".h5ad" )
    perturbed_expression_data_train_i, perturbed_expression_data_heldout_i = get_current_data_split(i, verbose = True)
    if args.amount_to_do in {"models", "missing_models"}:
        # Fit models!!
        if \
            (args.amount_to_do in {"models"}) or \
            (args.amount_to_do in {"missing_models"} and not os.path.isfile(h5ad)):
            try:
                os.unlink(h5ad)
            except FileNotFoundError:
                pass
            try:
                print(f"Fitting model for condition {i}")
                print(conditions.loc[i,:].T)
                grn = experimenter.do_one_run(
                    conditions = conditions, 
                    i = i,
                    train_data = perturbed_expression_data_train_i, 
                    test_data  = perturbed_expression_data_heldout_i,
                    networks = networks, 
                    outputs = outputs,
                    metadata = metadata,
                    human_tfs = DEFAULT_HUMAN_TFs,
                )
            except Exception as e: 
                if args.skip_bad_runs:
                    print(f"Caught exception {repr(e)} on experiment {i}; skipping.")
                else:
                    raise e
                continue

            if args.save_models:
                print("Saving models...", flush = True)
                grn.save_models( models )
            
            # Make predictions on test and (maybe) train set
            print("Generating predictions...", flush = True)
            if conditions.loc[i, "starting_expression"] == "control":
                predictions       = None
                predictions_train = None
            elif conditions.loc[i, "starting_expression"] == "heldout":
                print("Setting up initial conditions.")
                predictions       = perturbed_expression_data_heldout_i.copy()
                predictions_train = perturbed_expression_data_train_i.copy()
            else:
                raise ValueError(f"Unexpected value of 'starting_expression' in metadata: { conditions.loc[i, 'starting_expression'] }")
            print("Running GRN.predict()...")
            predictions   = grn.predict(
                [
                    (r[1][0], r[1][1]) 
                    for r in perturbed_expression_data_heldout_i.obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
                ],
                predictions = predictions,
                control_subtype = conditions.loc[i, "control_subtype"], 
                feature_extraction_requires_raw_data = (grn.feature_extraction == "geneformer"),
            )
            predictions.obs.index = perturbed_expression_data_heldout_i.obs.index.copy()
            # Sometimes AnnData has trouble saving pandas bool columns and sets, and they aren't needed here anyway.
            try:
                del predictions.obs["is_control"] 
                del predictions.obs["is_treatment"] 
                predictions.uns["perturbed_and_measured_genes"] = list(predictions.uns["perturbed_and_measured_genes"])
                predictions.uns["perturbed_but_not_measured_genes"] = list(predictions.uns["perturbed_but_not_measured_genes"])
            except KeyError as e:
                pass
            print("Saving predictions...")
            experimenter.safe_save_adata( predictions, h5ad )
            del predictions
            if args.save_trainset_predictions:
                fitted_values = grn.predict(
                    [
                        (r[1][0], r[1][1]) 
                        for r in perturbed_expression_data_train_i.obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
                    ], 
                    predictions = predictions_train, 
                    feature_extraction_requires_raw_data = (grn.feature_extraction == "geneformer"),
                )
                fitted_values.obs.index = perturbed_expression_data_train_i.obs.index.copy()
                # Sometimes AnnData has trouble saving pandas bool columns, and they aren't needed here anyway.
                experimenter.safe_save_adata( fitted_values, h5ad_fitted )
            print("... done.", flush = True)
            del grn

# Evaluate the results
if args.amount_to_do in {"models", "missing_models", "evaluations"}:
    print("Retrieving saved predictions", flush = True)
    conditions = experimenter.load_successful_conditions(outputs)
    predictions = {i:sc.read_h5ad( os.path.join(outputs, "predictions",   str(i) + ".h5ad" ), backed='r' ) for i in conditions.index}
    try:
        fitted_values = {i:sc.read_h5ad( os.path.join(outputs, "fitted_values", str(i) + ".h5ad" ), backed='r' ) for i in conditions.index}
    except FileNotFoundError:
        fitted_values = None
    for i in conditions.index:
        perturbed_expression_data_train_i, perturbed_expression_data_heldout_i = get_current_data_split(i)
        try:
            assert predictions[i].shape[0]==perturbed_expression_data_heldout_i.shape[0]
        except AssertionError:
            print(f"Object shapes for condition {i}: (observed, predicted):", flush = True)
            print((predictions[i].shape, perturbed_expression_data_heldout_i.shape), flush = True)
            raise AssertionError("Predicted and observed anndata are different shapes.")
        assert all(
                    np.sort(predictions[i].obs_names) == np.sort(perturbed_expression_data_heldout_i.obs_names)
                ), f"For condition {i}, set of observations is different between observed and predicted."        

    print("(Re)doing evaluations")
    evaluationPerPert, evaluationPerTarget = evaluator.evaluateCausalModel(
        get_current_data_split = get_current_data_split, 
        predicted_expression =  predictions,
        is_test_set = True,
        conditions = conditions,
        outputs = outputs,
        classifier = None,
        do_scatterplots = False,
    )
    evaluationPerPert.to_parquet(   os.path.join(outputs, "evaluationPerPert.parquet"))
    evaluationPerTarget.to_parquet( os.path.join(outputs, "evaluationPerTarget.parquet"))
    if fitted_values is not None:
        print("(Re)doing evaluations on (training set predictions)")
        evaluationPerPertTrainset, evaluationPerTargetTrainset = evaluator.evaluateCausalModel(
            get_current_data_split = get_current_data_split, 
            predicted_expression =  fitted_values,
            is_test_set = False,
            conditions = conditions,
            outputs = os.path.join(outputs, "trainset_performance"),
            classifier = None,
            do_scatterplots = False,
        )
        os.makedirs(os.path.join(outputs, "trainset_performance"), exist_ok=True)
        evaluationPerPertTrainset.to_parquet(   os.path.join(outputs, "trainset_performance", "evaluationPerPert.parquet"))
        evaluationPerTargetTrainset.to_parquet( os.path.join(outputs, "trainset_performance", "evaluationPerTarget.parquet"))
        
# Plot the results
if args.amount_to_do in {"plots", "models", "missing_models", "evaluations"}:
    print("Retrieving saved predictions", flush = True)
    conditions = experimenter.load_successful_conditions(outputs)
    predictions   = {i:sc.read_h5ad( os.path.join(outputs, "predictions",   str(i) + ".h5ad" ) ) for i in conditions.index}
    try:
        fitted_values = {i:sc.read_h5ad( os.path.join(outputs, "fitted_values", str(i) + ".h5ad" ) ) for i in conditions.index}
    except FileNotFoundError:
        fitted_values = None

    print("Retrieving saved evaluations", flush = True)
    evaluationPerPert   = pd.read_parquet(os.path.join(outputs, "evaluationPerPert.parquet"))
    evaluationPerTarget = pd.read_parquet(os.path.join(outputs, "evaluationPerTarget.parquet"))

    print("Plotting main summaries of results.")
    evaluator.makeMainPlots(
        evaluationPerPert, 
        evaluationPerTarget, 
        outputs = outputs, 
        factor_varied = metadata["factor_varied"],
        color_by = metadata["color_by"],
        facet_by=metadata["facet_by"],
    )
    try:
        evaluationPerPertTrainset   = pd.read_parquet(os.path.join(outputs, "trainset_performance", "evaluationPerPert.parquet"))
        evaluationPerTargetTrainset = pd.read_parquet(os.path.join(outputs, "trainset_performance", "evaluationPerTarget.parquet"))
        evaluator.makeMainPlots(
            evaluationPerPertTrainset, 
            evaluationPerTargetTrainset, 
            outputs = os.path.join(outputs, "trainset_performance"), 
            factor_varied = metadata["factor_varied"],
            color_by = metadata["color_by"],
            facet_by = metadata["facet_by"],
        )
    except FileNotFoundError:
        pass

print("Experiment done at " + str(datetime.datetime.now()), flush = True)