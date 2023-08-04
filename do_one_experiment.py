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
parser.add_argument("--skip_bad_runs", help="If true, keep running when some runs hit errors.", default = True, type = bool)
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
# For interactive use
if args.experiment_name is None:
    args = Namespace(**{
        "experiment_name": "1.4.2_2",
        "amount_to_do": "missing_models",
        "save_trainset_predictions": True,
        "save_models": False,
        "skip_bad_runs": False,
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
os.makedirs(os.path.join( outputs, "predictions"   ), exist_ok=True) 
os.makedirs(os.path.join( outputs, "fitted_values" ), exist_ok=True) 

for i in experiments.index:
    models      = os.path.join( outputs, "models",        str(i) )
    h5ad        = os.path.join( outputs, "predictions",   str(i) + ".h5ad" )
    h5ad_fitted = os.path.join( outputs, "fitted_values", str(i) + ".h5ad" )
    perturbed_expression_data = experimenter.filter_genes(perturbed_expression_data, num_genes = experiments.loc[i, "num_genes"], outputs = outputs)
    perturbed_expression_data_train[i], perturbed_expression_data_heldout[i] = experimenter.splitDataWrapper(
        perturbed_expression_data,
        networks = networks, 
        desired_heldout_fraction = experiments.loc[i, "desired_heldout_fraction"],  
        type_of_split            = experiments.loc[i, "type_of_split"],
        data_split_seed          = experiments.loc[i, "data_split_seed"],
    )
    
    # Delete the new copies and replace with shallow copies unless the data split is actually different
    def is_equal_anndata(ad1, ad2):
        return np.array_equal(ad1.X, ad2.X) and ad1.obs.equals(ad2.obs) and ad1.var.equals(ad2.var)
    if i>0 and \
    is_equal_anndata(perturbed_expression_data_train[i], perturbed_expression_data_train[i-1]) and \
    is_equal_anndata(perturbed_expression_data_heldout[i], perturbed_expression_data_heldout[i-1]):
        perturbed_expression_data_train[i] = perturbed_expression_data_train[i-1]
        perturbed_expression_data_heldout[i] = perturbed_expression_data_heldout[i-1]
        gc.collect()

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
                grn = experimenter.do_one_run(
                    experiments = experiments, 
                    i = i,
                    train_data = perturbed_expression_data_train[i], 
                    test_data  = perturbed_expression_data_heldout[i],
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
            if experiments.loc[i, "starting_expression"] == "control":
                starting_expression       = None
                starting_expression_train = None
            elif experiments.loc[i, "starting_expression"] == "heldout":
                starting_expression       = perturbed_expression_data_heldout[i].copy()
                starting_expression_train = perturbed_expression_data_train[i].copy()
            else:
                raise ValueError(f"Unexpected value of 'starting_expression' in metadata: { experiments.loc[i, 'starting_expression'] }")
            predictions   = grn.predict(
                [
                    (r[1][0], r[1][1]) 
                    for r in perturbed_expression_data_heldout[i].obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
                ], 
                starting_expression = starting_expression,
                control_subtype = experiments.loc[i, "control_subtype"]
            )   
            predictions.obs.index = perturbed_expression_data_heldout[i].obs.index.copy()
            # Sometimes AnnData has trouble saving pandas bool columns and sets, and they aren't needed here anyway.
            try:
                del predictions.obs["is_control"] 
                del predictions.obs["is_treatment"] 
                predictions.uns["perturbed_and_measured_genes"] = list(predictions.uns["perturbed_and_measured_genes"])
                predictions.uns["perturbed_but_not_measured_genes"] = list(predictions.uns["perturbed_but_not_measured_genes"])
            except KeyError as e:
                pass
            experimenter.safe_save_adata( predictions, h5ad )
            if args.save_trainset_predictions:
                fitted_values = grn.predict(
                    [
                        (r[1][0], r[1][1]) 
                        for r in perturbed_expression_data_train[i].obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
                    ], 
                    starting_expression = starting_expression_train
                )
                fitted_values.obs.index = perturbed_expression_data_train[i].obs.index.copy()
                # Sometimes AnnData has trouble saving pandas bool columns, and they aren't needed here anyway.
                experimenter.safe_save_adata( fitted_values, h5ad_fitted )
            print("... done.", flush = True)
            del grn

if args.amount_to_do in {"models", "missing_models", "evaluations"}:
    print("Retrieving saved predictions", flush = True)
    experiments = experimenter.load_successful_experiments(outputs)
    predictions = {i:sc.read_h5ad( os.path.join(outputs, "predictions",   str(i) + ".h5ad" ), backed='r' ) for i in experiments.index}
    try:
        fitted_values = {i:sc.read_h5ad( os.path.join(outputs, "fitted_values", str(i) + ".h5ad" ), backed='r' ) for i in experiments.index}
    except FileNotFoundError:
        fitted_values = None
    op = zip(predictions.values(), perturbed_expression_data_heldout.values())
    try:
        assert all(
            [obs.shape[0] == pred.shape[0] 
            for obs,pred in op]
        )
    except AssertionError:
        print("Object shapes (observed, predicted):", flush = True)
        op = zip(perturbed_expression_data_heldout.values(),predictions.values())
        print([(obs.shape, pred.shape) for obs,pred in op], flush = True)
        raise AssertionError("Predicted and observed anndata are different shapes.")
    assert all(
        [
            all(
                np.sort(predictions[i].obs_names) == np.sort(perturbed_expression_data_heldout[i].obs_names)
            ) 
            for i in predictions.keys()
        ]
    ) 
    
    print("(Re)doing evaluations")
    evaluationPerPert, evaluationPerTarget = evaluator.evaluateCausalModel(
        heldout = perturbed_expression_data_heldout,
        predictions = predictions,
        baseline = {
            i: perturbed_expression_data_train[i][[bool(b) for b in perturbed_expression_data_train[i].obs["is_control"]], :] 
            for i in perturbed_expression_data_train.keys()
        },
        experiments = experiments,
        outputs = outputs,
        classifier = None,
        do_scatterplots = False,
    )
    evaluationPerPert.to_parquet(   os.path.join(outputs, "evaluationPerPert.parquet"))
    evaluationPerTarget.to_parquet( os.path.join(outputs, "evaluationPerTarget.parquet"))
    if fitted_values is not None:
        print("(Re)doing evaluations on (training set predictions)")
        evaluationPerPertTrainset, evaluationPerTargetTrainset = evaluator.evaluateCausalModel(
            heldout = perturbed_expression_data_train,
            predictions = fitted_values,
            baseline = {
                i: perturbed_expression_data_train[i][perturbed_expression_data_train[i].obs["is_control"].astype(bool), :] 
                for i in perturbed_expression_data_train.keys()
            },
            experiments = experiments,
            outputs = os.path.join(outputs, "trainset_performance"),
            classifier = None,
            do_scatterplots = False,
        )
        os.makedirs(os.path.join(outputs, "trainset_performance"), exist_ok=True)
        evaluationPerPertTrainset.to_parquet(   os.path.join(outputs, "trainset_performance", "evaluationPerPert.parquet"))
        evaluationPerTargetTrainset.to_parquet( os.path.join(outputs, "trainset_performance", "evaluationPerTarget.parquet"))
        

if args.amount_to_do in {"plots", "models", "missing_models", "evaluations"}:
    print("Retrieving saved predictions", flush = True)
    experiments = experimenter.load_successful_experiments(outputs)
    predictions   = {i:sc.read_h5ad( os.path.join(outputs, "predictions",   str(i) + ".h5ad" ) ) for i in experiments.index}
    try:
        fitted_values = {i:sc.read_h5ad( os.path.join(outputs, "fitted_values", str(i) + ".h5ad" ) ) for i in experiments.index}
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

    print("Studying predictability for each target gene.")
    evaluationPerTarget = evaluator.studyPredictableGenes(
        evaluationPerTarget, 
        train_data = perturbed_expression_data_train[0], 
        test_data = perturbed_expression_data_heldout[0], 
        save_path = outputs,
        factor_varied = metadata["factor_varied"],
        genes_considered_as = "targets"
    )
    evaluationPerTarget.to_parquet(f"{outputs}/evaluationPerTarget.parquet")
    print("Studying predictability for each perturbation.")
    evaluationPerPert = evaluator.studyPredictableGenes(
        evaluationPerTarget = evaluationPerPert, 
        train_data = perturbed_expression_data_train[0], 
        test_data = perturbed_expression_data_heldout[0], 
        save_path = outputs,
        factor_varied = metadata["factor_varied"],        
        genes_considered_as = "perturbations"
    )
    evaluationPerPert.to_parquet(f"{outputs}/evaluationPerPert.parquet")

    print("Plotting all data and predictions for some example target genes.")
    if fitted_values is not None:
        for type in ["best", "worst", "random"]:
            if type=="best":
                genes = evaluationPerTarget.nlargest( 5, "mae_benefit")["target"]
            elif type=="worst":
                genes = evaluationPerTarget.nsmallest(5, "mae_benefit")["target"]
            else:
                genes = np.random.choice(size=5, a=evaluationPerTarget["target"].unique())
            for gene in genes:
                print(gene)
                evaluator.plotOneTargetGene(
                    gene=gene, 
                    outputs=os.path.join(outputs,"targets", type), 
                    experiments=experiments,
                    factor_varied=metadata["factor_varied"],
                    train_data=perturbed_expression_data_train, 
                    heldout_data=perturbed_expression_data_heldout, 
                    fitted_values=fitted_values, 
                    predictions=predictions)
    


print("Experiment done at " + str(datetime.datetime.now()), flush = True)