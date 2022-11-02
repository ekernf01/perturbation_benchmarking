"""evaluator.py is a collection of functions for making and testing predictions about expression fold change after genetic perturbations.
It dedicates particular attention to interfacing with CellOracle, a thought-provoking and flexible perturbation prediction method.
"""
from ast import If
from multiprocessing.sharedctypes import Value
import celloracle as co
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gc
import anndata
from scipy.stats import spearmanr as spearmanr

import sys
import importlib
sys.path.append("src")
import predict
importlib.reload(predict)

def makeMainPlots(evaluationResults, factor_varied:str, outputs: str, default_level=None):
    """Redo the main plots summarizing an experiment.
    Args:
        factor_varied (str): Plots are automatically stratified based on this column of "experiments". 
        evaluationResults (pd.DataFrame, optional): By default, looks for results saved as a csv.
        outputs (str): put plots here
        default_level: a value found in evaluationResults.loc[:,factor_varied] to compare all others against. This will be used to remove the variation from gene to gene. 
    """
    stripchartMainFig = sns.stripplot(y = "spearman", 
                                                x=factor_varied, 
                                                data = evaluationResults)
    stripchartMainFig.set(ylabel="Spearman correlation \npredicted log fold change vs observed")
    stripchartMainFig.set_xticklabels(stripchartMainFig.get_xticklabels(), rotation=90)
    stripchartMainFig.figure.savefig(f'{outputs}/stripchart.pdf', bbox_inches="tight")
    plt.show()
    plt.figure()
    meanSEPlot = {}
    for readout in "spearman", "cell_fate_correct":
        meanSEPlot[readout] = sns.pointplot(factor_varied, y=readout, data=evaluationResults, dodge=True, join=False)
        meanSEPlot[readout].set(title="Mean and SE of " + readout)
        meanSEPlot[readout].set_xticklabels(meanSEPlot[readout].get_xticklabels(), rotation=90)
        meanSEPlot[readout].figure.savefig(f'{outputs}/MeanSEPlot{readout}.pdf', bbox_inches="tight")
        plt.show()
        plt.figure()
    return stripchartMainFig, meanSEPlot


def evaluateCausalModel(
    heldout:anndata.AnnData, 
    predictions:anndata.AnnData, 
    baseline:anndata.AnnData,
    experiments: pd.DataFrame, 
    outputs: str, 
    factor_varied: str,
    default_level, 
    classifier = None):
    """Compile plots and tables comparing heldout data and predictions for same. 

    Args:
        heldout (AnnData): Test data in the format specified by this project's collection of perturbation data.
        predictions: dictionary with keys equal to index of experiments. 
            Each value is an AnnData object.
        baseline (AnnData): Expression before perturbation, for use in calculating log fold change. 
        classifier (sklearn.LogisticRegression): Optional, to judge results on cell type accuracy. 
        experiments (pd.DataFrame): Metadata for the different combinations used in this experiment. 
        factor_varied (String): Plots are automatically stratified based on this column of "experiments". 
        default_level: a value found in experiments.loc[:,factor_varied] to compare all others against. This is for use in repeated-measures analyses to remove the variation from gene to gene. 
        outputs (String): Saves output here.
    """
    # Get spearman and classifier accuracy 
    evaluationResults = {}
    shared_var_names = list(set.intersection(*[set(predictions[experiment].var_names) for experiment in predictions.keys()]))
    for experiment in predictions.keys(): 
        evaluationResults[experiment] = \
            evaluateOnePrediction(
                expression =                        heldout[:, shared_var_names], 
                predictedExpression=predictions[experiment][:, shared_var_names],   
                baseline =                         baseline[:, shared_var_names],     
                doPlots=False, 
                classifier=classifier
            )[0]
        evaluationResults[experiment]["index"] = experiment
    evaluationResults = pd.concat(evaluationResults)
    evaluationResults = evaluationResults.merge(experiments, how = "left", right_index = True, left_on = "index")
    evaluationResults = pd.DataFrame(evaluationResults.to_dict())
    # TO DO: Estimate a per-gene fixed effect 
    # evaluationResults["spearman_gene_effect"] = 
    # Mark anything where predictions were unavailable
    noPredictionMade = evaluationResults.iloc[[x==0 for x in evaluationResults["spearman"]],:]['perturbation']
    noPredictionMade = set(noPredictionMade)
    noPredictionMade
    evaluationResults["somePredictionRefused"] = evaluationResults["perturbation"].isin(noPredictionMade) 
    # Save main results
    evaluationResults.to_csv(outputs +"/evaluationResults.csv")
    stripchartMainFig, meanSEPlot = makeMainPlots(factor_varied=factor_varied,  outputs=outputs, evaluationResults=evaluationResults)
    # Which genes are best/worst?
    hardest = evaluationResults.loc[evaluationResults["spearman"].idxmax(),"perturbation"]
    easiest = evaluationResults.loc[evaluationResults["spearman"].idxmin(),"perturbation"]
    evaluationResults.loc[evaluationResults["perturbation"]==hardest,:].to_csv(outputs +"/hardest.csv")
    evaluationResults.loc[evaluationResults["perturbation"]==easiest,:].to_csv(outputs +"/easiest.csv")
    return evaluationResults, stripchartMainFig, meanSEPlot

def evaluateOnePrediction(expression: anndata.AnnData, predictedExpression: anndata.AnnData, baseline: anndata.AnnData, doPlots=False, classifier = None, do_careful_checks = True):
    '''Compare observed against predicted, for expression, fold-change, or cell type.

            Parameters:
                    expression (AnnData): 
                        the observed expression post-perturbation (log-scale in expression.X). 
                    predictedExpression (AnnData): 
                        the cellOracle prediction (log-scale). Elements of predictedExpression.X may be np.nan for 
                        missing predictions, often one gene missing from all samples or one sample missing for all genes.
                        predictedExpression.obs must contain columns "perturbation" (symbol of targeted gene) 
                        and "expression_level_after_perturbation" (e.g. 0 for knockouts). 
                    baseline (AnnData): 
                        control expression level (log-scale)
                    classifier (sklearn logistic regression classifier): 
                        optional machine learning classifier to assign cell fate. 
                        Must have a predict() method capable of taking a value from expression or predictedExpression and returning a single class label. 
                    do_careful_checks (bool): check gene name and expression level associated with each perturbation.
                        They must match between expression and predictionExpression.
            Returns:
                    Pandas DataFrame with Spearman correlation between predicted and observed 
                    log fold change over control.
    '''
    "log fold change using Spearman correlation and (optionally) cell fate classification."""
    if not expression.X.shape == predictedExpression.X.shape:
        raise ValueError("expression and predictedExpression must have the same shape.")
    if not expression.X.shape[1] == baseline.X.shape[1]:
        raise ValueError("expression and baseline must have the same number of genes.")
    predictedExpression.obs_names = expression.obs_names
    baseline = baseline.X.mean(axis=0).squeeze()
    plots = {}
    metrics = pd.DataFrame(index = predictedExpression.obs.index, columns = ["spearman", "spearmanp", "cell_fate_correct"])
    for pert in predictedExpression.obs.index:
        if do_careful_checks:
            assert all(
                                 expression.obs.loc[pert, ["perturbation", "expression_level_after_perturbation"]] == \
                        predictedExpression.obs.loc[pert, ["perturbation", "expression_level_after_perturbation"]] 
                    )
        observed  = expression[         pert,:].X.squeeze()
        predicted = predictedExpression[pert,:].X.squeeze()
        if type(predicted) is float and np.isnan(predicted):
            metrics.loc[pert,["spearman","spearmanp", "cell_fate_correct"]] = 0,1,np.nan
        else:
            metrics.loc[pert,["spearman","spearmanp"]] = spearmanr(observed - baseline, predicted - baseline)
            if classifier is not None:
                class_observed  = classifier.predict(np.reshape(observed,  (1, -1)))[0]
                class_predicted = classifier.predict(np.reshape(predicted, (1, -1)))[0]
                metrics.loc[pert,"cell_fate_correct"] = 1.0*(class_observed==class_predicted)            
        if doPlots:
            plots[pert] = sns.scatterplot(x=observed, y=predicted)
            plots[pert].set(title=pert + " (Spearman rho="+ str(round(metrics.loc[pert,"spearman"])) +")")
            plots[pert].set_xlabel("Observed log fc", fontsize = 20)
            plots[pert].set_ylabel("Predicted log fc", fontsize = 20)
            plots[pert].plot()
            plt.figure()
    metrics["perturbation"] = metrics.index
    return metrics, plots                 
    
def networkEdgesToMatrix(networkEdges, regulatorColumn=0, targetColumn=1):
    """Reformat a network from a two-column dataframe to the way that celloracle needs its input."""
    X = pd.crosstab(networkEdges.iloc[:,targetColumn], networkEdges.iloc[:,regulatorColumn])
    del networkEdges
    gc.collect()
    X = 1.0*(X > 0)
    X = X.rename_axis('gene_short_name').reset_index()
    X = X.rename_axis('peak_id').reset_index()
    X = makeNetworkSparse(X, 0.0)
    gc.collect()
    return X


def countMatrixEdges(network):
    """Count the number of connections in a network that is formatted how CO expects it to be formatted."""
    return 1.0*network.iloc[:,2:].sum().sum()

humanTFs = pd.read_csv("../accessory_data/humanTFs.csv")
targetGenes = co.data.load_human_promoter_base_GRN()["gene_short_name"]

def pivotNetworkWideToLong(network_wide: pd.DataFrame):
    """Convert from CellOracle's preferred format to a triplet format

    Args:
        network_wide (pd.DataFrame): GRN structure in CellOracle's usual format
    """
    network_long = pd.concat([
        pd.DataFrame({
            "regulator": tf,
            "target": network_wide.loc[network_wide[tf]==1, "gene_short_name"],
            "weight": 1,
        })
        for tf in network_wide.columns[2:]
    ])
    return network_long

def makeRandomNetwork(density = 0, seed = 0, TFs = humanTFs['HGNC symbol'], targetGenes = targetGenes ):
    """Generate a random network formatted the way that celloracle needs its input."""
    np.random.seed(seed)
    X = pd.DataFrame(
            np.random.binomial(
                n = 1, 
                p=density,
                size=(
                    len(targetGenes), 
                    len(TFs)
                )
            ),
            columns = TFs, 
            index = targetGenes
        )
    X.rename_axis('gene_short_name', inplace=True)
    X.reset_index(inplace=True)
    X.rename_axis('peak_id', inplace=True)
    X.reset_index(inplace=True)
    # CellOracle's preferred format wastes gobs of memory unless you sparsify.
    X = makeNetworkSparse(X, round(density))
    gc.collect()
    return X


def makeNetworkSparse(X, defaultValue):
    """Save memory by making a sparse representation of a base network"""
    X.iloc[:,2:] = X.iloc[:,2:].astype(pd.SparseDtype("float", defaultValue))
    return X


def makeNetworkDense(X):
    """Undo makeNetworkSparse"""
    X.iloc[:, 2:] = np.array(X.iloc[:, 2:])   #undo sparse representation         
    return X


def splitData(adata, allowedRegulators, minTestSetSize = 250):
    """Determine a train-test split satisfying constraints imposed by base networks and available data.
    
A few factors complicate the training-test split. 

- Perturbed genes may be absent from most base GRN's due to lack of motif information or ChIP data. These are excluded from the test data to avoid obvious failure cases.
- Perturbed genes may not be measured. These are excluded from test data because we don't know to what extent they were overexpressed.

In both cases, we still use those perturbed profiles as training data, hoping they will provide useful info about attainable cell states and downstream causal effects. 

For some collections of base networks, there are many factors ineligible for use as test data -- so many that we use all the eligible ones for test and the only ineligible ones for training. For other cases, such as dense base networks, we have more flexibility, so we send some perturbations to the training set at random even if we would be able to use them in the test set.

parameters:
    - adata: AnnData object satisfying the expectations outlined in the accompanying collection of perturbation data.
    - allowedRegulators: list or set of features allowed to be in the test set. In CellOracle, this is usually constrained by motif/chip availability. 

    """
    allowedRegulators = set(allowedRegulators)
    allowedRegulators = allowedRegulators.intersection(adata.uns["perturbed_and_measured_genes"])
    if len(allowedRegulators) <= minTestSetSize:
        raise ValueError(f"minTestSetSize was set to {minTestSetSize} but only {len(allowedRegulators)} perturbed conditions are available.")
    testSetPerturbations     = set(adata.obs["perturbation"]).intersection(allowedRegulators)
    trainingSetPerturbations = set(adata.obs["perturbation"]).difference(allowedRegulators)
    if len(trainingSetPerturbations) < minTestSetSize:
        swap = np.random.default_rng(seed=0).choice(list(testSetPerturbations), 
                                               minTestSetSize - len(trainingSetPerturbations), 
                                               replace = False)
        testSetPerturbations = testSetPerturbations.difference(swap)
        trainingSetPerturbations = trainingSetPerturbations.union(swap)
    adata_train    = adata[adata.obs["perturbation"].isin(trainingSetPerturbations),:]
    adata_heldout  = adata[adata.obs["perturbation"].isin(testSetPerturbations),    :]
    adata_train.uns[  "perturbed_and_measured_genes"] = adata_train.uns[  "perturbed_and_measured_genes"].intersection(trainingSetPerturbations)
    adata_heldout.uns["perturbed_and_measured_genes"] = adata_heldout.uns["perturbed_and_measured_genes"].intersection(testSetPerturbations)
    adata_train.uns[  "perturbed_but_not_measured_genes"] = adata_train.uns[  "perturbed_but_not_measured_genes"].intersection(trainingSetPerturbations)
    adata_heldout.uns["perturbed_but_not_measured_genes"] = adata_heldout.uns["perturbed_but_not_measured_genes"].intersection(testSetPerturbations)
    print("Test set size:")
    print(len(testSetPerturbations))
    print("Training set size:")
    print(len(trainingSetPerturbations))
    return adata_train, adata_heldout


def downsample(adata: anndata.AnnData, proportion: float, seed = None, proportion_genes = 1):
    """Downsample training data to a given fraction, always keeping controls. 
    Args:
        adata (anndata.AnnData): _description_
        proportion (float): fraction of observations to keep. You may end up with a little extra because all controls are kept.
        proportion_genes (float): fraction of cells to keep. You may end up with a little extra because all perturbed genes are kept.
        seed (_type_, optional): RNG seed. Seed defaults to proportion so if you ask for 80% of cells, you get the same 80% every time.

    Returns:
        anndata.AnnData: Subsampled data.
    """
    if seed is None:
        seed = proportion
    np.random.seed(int(np.round(seed)))
    mask       = np.random.choice(a=[True, False], size=adata.obs.shape[0], p=[proportion,       1-proportion], replace = True)
    mask_genes = np.random.choice(a=[True, False], size=adata.var.shape[0], p=[proportion_genes, 1-proportion_genes], replace = True)
    adata = adata[adata.obs["is_control"] | mask, :].copy()
    perturbed_genes_remaining = set(adata.obs["perturbation"])
    adata = adata[:, [adata.var.index.isin(perturbed_genes_remaining)] | mask_genes].copy()
    print(adata.obs.shape)
    adata.uns["perturbed_but_not_measured_genes"] = set(adata.obs["perturbation"]).difference(  set(adata.var_names))
    adata.uns["perturbed_and_measured_genes"]     = set(adata.obs["perturbation"]).intersection(set(adata.var_names))
    return adata