import os
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
import pandas as pd
import scanpy as sc
import gc
import argparse
from argparse import Namespace
import datetime
# Deal with various data and modules specific to this project
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
import evaluator
import experimenter
import load_perturbations
import ggrn
importlib.reload(evaluator)
importlib.reload(experimenter)
importlib.reload(load_perturbations)
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbation_data/perturbations"
# User input: name of experiment and whether to fully rerun or just remake plots. 
parser = argparse.ArgumentParser("experimenter")
parser.add_argument("--experiment_name", help="Unique id for the experiment.", type=str)
parser.add_argument("--save_models",     help="If true, save model objects.", default = False, action = "store_true")
parser.add_argument("--save_trainset_predictions", help="If true, make & save predictions of training data.", default = False, action = "store_true")
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
        "experiment_name":"1.0_0",
        "amount_to_do": "missing_models",
        "save_trainset_predictions": True,
        "save_models": False,
    })
# Additional inputs and outputs
print("Running experiment", flush = True)
outputs = os.path.join("experiments", args.experiment_name, "outputs")
os.makedirs(outputs, exist_ok=True)
metadata = experimenter.validate_metadata(experiment_name=args.experiment_name)


print("Starting at " + str(datetime.datetime.now()), flush = True)
perturbed_expression_data, networks, experiments = experimenter.set_up_data_networks_conditions(
    metadata,
    amount_to_do = args.amount_to_do, 
    outputs = outputs,
)

perturbed_expression_data_train = {}
perturbed_expression_data_heldout = {}
if args.amount_to_do in {"models", "missing_models", "evaluations"}:
    os.makedirs(os.path.join( outputs, "predictions"   ), exist_ok=True) 
    os.makedirs(os.path.join( outputs, "fitted_values" ), exist_ok=True) 
    for i in experiments.index:
        models      = os.path.join( outputs, "models",        str(i) )
        h5ad        = os.path.join( outputs, "predictions",   str(i) + ".h5ad" )
        h5ad_fitted = os.path.join( outputs, "fitted_values", str(i) + ".h5ad" )
        perturbed_expression_data = experimenter.filter_genes(perturbed_expression_data, num_genes = experiments.loc[i, "num_genes"])
        perturbed_expression_data_train[i], perturbed_expression_data_heldout[i] = experimenter.splitDataWrapper(
            perturbed_expression_data,
            networks = networks, 
            desired_heldout_fraction = experiments.loc[i, "desired_heldout_fraction"],  
            type_of_split            = experiments.loc[i, "type_of_split"],
            data_split_seed          = experiments.loc[i, "data_split_seed"],
        )
        if \
            (args.amount_to_do in {"models"}) or \
            (args.amount_to_do in {"missing_models"} and not os.path.isfile(h5ad)):
            
            try:
                os.unlink(h5ad)
            except FileNotFoundError:
                pass
            try:
                grn = experimenter.do_one_run(
                    experiments, 
                    i,
                    train_data = perturbed_expression_data_train[i], 
                    test_data  = perturbed_expression_data_heldout[i],
                    networks = networks, 
                    outputs = outputs,
                    metadata = metadata,
                )
            except Exception as e:
                print(f"Caught exception {repr(e)} on experiment {i}; skipping.")
                continue

            if args.save_models:
                print("Saving models...", flush = True)
                grn.save_models( models )
            print("Saving predictions...", flush = True)
            predictions   = grn.predict([
                (r[1][0], r[1][1]) 
                for r in perturbed_expression_data_heldout[i].obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
            ])   
            predictions.obs.index = perturbed_expression_data_heldout[i].obs.index.copy()
            predictions.write_h5ad( h5ad )
            if args.save_trainset_predictions:
                fitted_values = grn.predict([
                    (r[1][0], r[1][1]) 
                    for r in perturbed_expression_data_train[i].obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
                ])
                fitted_values.obs.index = perturbed_expression_data_train[i].obs.index.copy()
                fitted_values.write_h5ad( h5ad_fitted )
            print("... done.", flush = True)
            del grn

if args.amount_to_do in {"models", "missing_models", "evaluations"}:
    print("Retrieving saved predictions", flush = True)
    experiments = pd.read_csv( os.path.join(outputs, "experiments.csv") )
    predictions = {i:sc.read_h5ad( os.path.join(outputs, "predictions",   str(i) + ".h5ad" ) ) for i in experiments.index}
    try:
        fitted_values = {i:sc.read_h5ad( os.path.join(outputs, "fitted_values", str(i) + ".h5ad" ) ) for i in experiments.index}
    except FileNotFoundError:
        fitted_values = None
    assert all(
        [obs.shape[0] == pred.shape[0] 
        for obs,pred in zip(predictions.values(), perturbed_expression_data_heldout.values())]
    )
    assert all(
        [all(
            np.sort(predictions[i].obs_names) == np.sort(perturbed_expression_data_heldout[i].obs_names)) 
            for i in range(len(predictions.keys())
        )]
    )
    np.sort(predictions[0].obs_names)
    np.sort(perturbed_expression_data_heldout[0].obs_names)

    evaluationPerPert, evaluationPerTarget = evaluator.evaluateCausalModel(
        heldout = perturbed_expression_data_heldout,
        predictions = predictions,
        baseline = {
            i: perturbed_expression_data_train[i][perturbed_expression_data_train[i].obs["is_control"], :] 
            for i in perturbed_expression_data_train.keys()
        },
        experiments = experiments,
        outputs = outputs,
        factor_varied = metadata["factor_varied"],
        default_level = None,
        classifier = None,
        do_scatterplots = False,
    )
    evaluationPerPert.to_parquet(   os.path.join(outputs, "evaluationPerPert.parquet"))
    evaluationPerTarget.to_parquet( os.path.join(outputs, "evaluationPerTarget.parquet"))
    if fitted_values is not None:
        evaluationPerPertTrainset, evaluationPerTargetTrainset = evaluator.evaluateCausalModel(
            heldout = perturbed_expression_data_train,
            predictions = fitted_values,
            baseline = {
                i: perturbed_expression_data_train[i][perturbed_expression_data_train[i].obs["is_control"], :] 
                for i in perturbed_expression_data_train.keys()
            },
            experiments = experiments,
            outputs = os.path.join(outputs, "trainset_performance"),
            factor_varied = metadata["factor_varied"],
            default_level = None,
            classifier = None,
            do_scatterplots = False,
        )
        evaluationPerPertTrainset.to_parquet(   os.path.join(outputs, "trainset_performance", "evaluationPerPert.parquet"))
        evaluationPerTargetTrainset.to_parquet( os.path.join(outputs, "trainset_performance", "evaluationPerTarget.parquet"))
        

if args.amount_to_do in {"plots", "models", "missing_models", "evaluations"}:
    print("Retrieving saved predictions", flush = True)
    experiments = pd.read_csv( os.path.join(outputs, "experiments.csv") )
    predictions   = {i:sc.read_h5ad( os.path.join(outputs, "predictions",   str(i) + ".h5ad" ) ) for i in experiments.index}
    try:
        fitted_values = {i:sc.read_h5ad( os.path.join(outputs, "fitted_values", str(i) + ".h5ad" ) ) for i in experiments.index}
    except FileNotFoundError:
        fitted_values = None
    print("Retrieving saved evaluations", flush = True)
    evaluationPerPert   = pd.read_parquet(os.path.join(outputs, "evaluationPerPert.parquet"))
    evaluationPerTarget = pd.read_parquet(os.path.join(outputs, "evaluationPerTarget.parquet"))
    evaluationPerTarget = evaluator.studyPredictableGenes(
        evaluationPerTarget, 
        perturbed_expression_data_train[0], 
        save_path = outputs,
        factor_varied = metadata["factor_varied"],
    )
    if fitted_values is not None:
        for type in ["best", "worst", "random"]:
            if type=="best":
                genes = evaluationPerTarget.nlargest( 100, "mae_benefit")["target"]
            elif type=="worst":
                genes = evaluationPerTarget.nsmallest(100, "mae_benefit")["target"]
            else:
                genes = np.random.choice(size=10, a=evaluationPerTarget["target"].unique())
            for gene in genes:
                evaluator.plotOneTargetGene(
                    gene=gene, 
                    outputs=os.path.join(outputs,"target_genes_best_worst", type), 
                    experiments=experiments,
                    factor_varied=metadata["factor_varied"],
                    train_data=perturbed_expression_data_train, 
                    heldout_data=perturbed_expression_data_heldout, 
                    fitted_values=fitted_values, 
                    predictions=predictions)
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