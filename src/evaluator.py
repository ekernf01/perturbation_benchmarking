"""evaluator.py is a collection of functions for making and testing predictions about expression fold change after genetic perturbations.
It dedicates particular attention to interfacing with CellOracle, a thought-provoking and flexible perturbation prediction method.
"""
from ast import If
from multiprocessing.sharedctypes import Value
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gc
import anndata
from scipy.stats import spearmanr as spearmanr
import os 
import shutil
import sys
import importlib
sys.path.append("src")
import altair as alt

def makeMainPlots(
    evaluationPerPert: pd.DataFrame, 
    evaluationPerTarget: pd.DataFrame, 
    outputs: str, 
    factor_varied:str, 
    facet_by: str = None, 
    color_by: str = None, 
    metrics = ['spearman', 'mse']
    ):
    """Redo the main plots summarizing an experiment.
    Args:
        evaluationPerPert (pd.DataFrame)
        evaluationPerTarget (pd.DataFrame)
        factor_varied (str): Plots are automatically colored based on this column of "evaluationPerPert". 
        facet_by (str): Plots are automatically stratified based on this column of "evaluationPerPert". 
        outputs (str): folder to save plots in
        default_level: a value found in evaluationPerPert.loc[:,factor_varied] to compare all others against. This will be used to remove the variation from gene to gene. 
        metrics: How to measure performance. 
    """
    # Sometimes the index is too complex for Altair to handle correctly (tuples)
    try:
        evaluationPerPert.index = [p[1] for p in evaluationPerPert.index]
    except IndexError:
        pass
    vlnplot = {}
    for metric in metrics:
        group_mean_by = [factor_varied]
        if facet_by is not None:
            group_mean_by.append(facet_by)
        if color_by is not None:
            group_mean_by.append(color_by)
        means = evaluationPerPert.groupby(group_mean_by, as_index=False)[[metric]].mean()
        vlnplot[metric] = alt.Chart(
            data = evaluationPerPert, 
            title = f"{metric} (predicted log fold change vs observed)"
        ).mark_boxplot(extent='min-max') + alt.Chart(data = means).mark_point(color="black")
        if color_by is not None:
            vlnplot[metric]=vlnplot[metric].encode(
                y=alt.Y(f'{metric}:Q'),
                color=color_by,
                x=alt.X(
                    factor_varied + ':N'
                )
            ).properties(
                width=400,
                height=400
            )
        else:
            vlnplot[metric] = vlnplot[metric].encode(
                y=alt.Y(f'{metric}:Q'),
                x=alt.X(
                    factor_varied + ':N'
                )
            ).properties(
                width=400,
                height=400
            )
        if facet_by is not None:
            vlnplot[metric] = vlnplot[metric].facet(
                facet_by + ':N',
                columns=int(np.ceil(np.sqrt(len(evaluationPerPert[facet_by].unique()))))
            )
        vlnplot[metric].save(f'{outputs}/{metric}.pdf')
    return vlnplot

def evaluateCausalModel(
    heldout:anndata.AnnData, 
    predictions:anndata.AnnData, 
    baseline:anndata.AnnData,
    experiments: pd.DataFrame, 
    outputs: str, 
    factor_varied: str,
    default_level, 
    classifier = None, 
    do_scatterplots = True):
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
    evaluationPerPert = {}
    evaluationPerTarget = {}
    shared_var_names = list(set.intersection(*[set(predictions[experiment].var_names) for experiment in predictions.keys()]))
    for experiment in predictions.keys(): 
        evaluationPerPert[experiment], evaluationPerTarget[experiment] = \
            evaluateOnePrediction(
                expression =                        heldout[:, shared_var_names], 
                predictedExpression=predictions[experiment][:, shared_var_names],   
                baseline =                         baseline[:, shared_var_names],     
                doPlots=do_scatterplots, 
                outputs = os.path.join(outputs, "plots", str(experiment)),
                classifier=classifier
            )
        evaluationPerPert[experiment]["index"]   = experiment
        evaluationPerTarget[experiment]["index"] = experiment
    evaluationPerPert   = pd.concat(evaluationPerPert)
    evaluationPerTarget = pd.concat(evaluationPerTarget)
    evaluationPerPert   = evaluationPerPert.merge(experiments,   how = "left", right_index = True, left_on = "index")
    evaluationPerTarget = evaluationPerTarget.merge(experiments, how = "left", right_index = True, left_on = "index")
    evaluationPerPert   = pd.DataFrame(evaluationPerPert.to_dict())
    evaluationPerTarget = pd.DataFrame(evaluationPerTarget.to_dict())
    # Mark anything where predictions were unavailable
    noPredictionMade = evaluationPerPert.iloc[[x==0 for x in evaluationPerPert["spearman"]],:]['perturbation']
    noPredictionMade = set(noPredictionMade)
    noPredictionMade
    evaluationPerPert["somePredictionRefused"] = evaluationPerPert["perturbation"].isin(noPredictionMade) 
    return evaluationPerPert, evaluationPerTarget

def evaluateOnePrediction(
    expression: anndata.AnnData, 
    predictedExpression: anndata.AnnData, 
    baseline: anndata.AnnData, 
    outputs,
    doPlots=False, 
    classifier = None, 
    do_careful_checks = True):
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
                    doPlots (bool): Make a scatterplot showing observed vs predicted, one dot per gene. 
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
    metrics = pd.DataFrame(index = predictedExpression.obs.index, columns = ["spearman", "spearmanp", "cell_fate_correct", "mse"])
    metrics_per_target = pd.DataFrame(index = predictedExpression.var.index, columns = ["mse"])
    for target in predictedExpression.var.index:
        observed  = expression[         :,target].X.squeeze()
        predicted = predictedExpression[:,target].X.squeeze()
        metrics_per_target.loc[target,["mse"]] = np.linalg.norm(observed - predicted)**2
    
    for pert in predictedExpression.obs.index:
        if do_careful_checks:
            assert all(
                                 expression.obs.loc[pert, ["perturbation", "expression_level_after_perturbation"]].fillna(0) == \
                        predictedExpression.obs.loc[pert, ["perturbation", "expression_level_after_perturbation"]].fillna(0) 
                    )
        observed  = expression[         pert,:].X.squeeze()
        predicted = predictedExpression[pert,:].X.squeeze()
        def is_constant(x):
            return np.std(x)<1e-12
        if type(predicted) is float and np.isnan(predicted) or is_constant(predicted - baseline) or is_constant(observed - baseline):
            metrics.loc[pert,["spearman","spearmanp", "cell_fate_correct", "mse"]] = 0,1,np.nan,np.nan
        else:
            metrics.loc[pert,["spearman","spearmanp"]] = spearmanr(observed - baseline, predicted - baseline)
            metrics.loc[pert,"mse"] = np.linalg.norm(observed - predicted)**2
            if classifier is not None:
                class_observed  = classifier.predict(np.reshape(observed,  (1, -1)))[0]
                class_predicted = classifier.predict(np.reshape(predicted, (1, -1)))[0]
                metrics.loc[pert,"cell_fate_correct"] = 1.0*(class_observed==class_predicted)  

    metrics["spearman"] = metrics["spearman"].astype(float)
    hardest = metrics["spearman"].idxmin()
    easiest = metrics["spearman"].idxmax()
    for pert in metrics.index:
        observed  = expression[         pert,:].X.squeeze()
        predicted = predictedExpression[pert,:].X.squeeze()
        is_hardest = hardest==pert
        is_easiest = easiest==pert
        if doPlots | is_hardest | is_easiest:
            os.makedirs(outputs, exist_ok = True)
            diagonal = alt.Chart(
                pd.DataFrame({
                    "x":[-1, 1],
                    "y":[-1,1 ],
                })
            ).mark_line(color= 'black').encode(
                x= 'x',
                y= 'y',
            )
            scatterplot = alt.Chart(
                pd.DataFrame({
                    "Observed log fc": observed-baseline, 
                    "Predicted log fc": predicted-baseline, 
                    "Baseline expression":baseline,
                })
            ).mark_circle().encode(
                x="Observed log fc:Q",
                y="Predicted log fc:Q",
                color="Baseline expression:Q",
            ).properties(
                title=pert + " (Spearman rho="+ str(round(metrics.loc[pert,"spearman"], ndigits=2)) +")"
            ) + diagonal
            scatterplot.save(os.path.join(outputs, f"{pert}.pdf"))
            if is_easiest:
                scatterplot.save(os.path.join(outputs, f"_easiest({pert}).pdf"))
            if is_hardest:
                scatterplot.save(os.path.join(outputs, f"_hardest({pert}).pdf"))
    metrics["perturbation"] = metrics.index
    return metrics, metrics_per_target
    
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

def makeRandomNetwork(targetGenes, density = 0, seed = 0, TFs = humanTFs['HGNC symbol'] ):
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


def splitData(adata, allowedRegulators, desired_heldout_fraction = 0.5):
    """Determine a train-test split satisfying constraints imposed by base networks and available data.
    
    A few factors complicate the training-test split. 

    - Perturbed genes may be absent from most base GRN's due to lack of motif information or ChIP data. These are 
        excluded from the test data to avoid obvious failure cases.
    - Perturbed genes may not be measured. These are excluded from test data because we don't know to what extent 
        they were overexpressed.

    In both cases, we still use those perturbed profiles as training data, hoping they will provide useful info about 
    attainable cell states and downstream causal effects. 

    For some collections of base networks, there are many factors ineligible for use as test data -- so many that 
    we use all the eligible ones for test and the only ineligible ones for training. 
    For other cases, such as dense base networks, we have more flexibility, so we send some perturbations to the 
    training set at random even if we would be able to use them in the test set.

    parameters:
    - adata: AnnData object satisfying the expectations outlined in the accompanying collection of perturbation data.
    - allowedRegulators: list or set of features allowed to be in the test set. In CellOracle, this is usually constrained 
        by motif/chip availability. 

    """
    # For a deterministic result when downsampling an iterable, setting a seed alone is not enough.
    # Must also avoid the use of sets. 
    get_unique_keep_order = lambda x: list(dict.fromkeys(x))
    allowedRegulators = [p for p in allowedRegulators if p in adata.uns["perturbed_and_measured_genes"]]
    testSetEligible   = [p for p in adata.obs["perturbation"] if p     in allowedRegulators]
    testSetIneligible = [p for p in adata.obs["perturbation"] if p not in allowedRegulators]
    allowedRegulators = get_unique_keep_order(allowedRegulators)
    testSetEligible   = get_unique_keep_order(testSetEligible)
    testSetIneligible = get_unique_keep_order(testSetIneligible)
    eligible_heldout_fraction = len(testSetEligible)/(0.0+len(allowedRegulators))
    if eligible_heldout_fraction < desired_heldout_fraction:
        print("Not enough profiles for the desired_heldout_fraction. Will use all available.")
        testSetPerturbations = testSetEligible
        trainingSetPerturbations = testSetIneligible
    elif eligible_heldout_fraction == desired_heldout_fraction: #nailed it
        testSetPerturbations = testSetEligible
        trainingSetPerturbations = testSetIneligible
    else:
        numExcessTestEligible = int(np.ceil((eligible_heldout_fraction - desired_heldout_fraction)*len(allowedRegulators)))
        excessTestEligible = np.random.default_rng(seed=0).choice(
            testSetEligible, 
            numExcessTestEligible, 
            replace = False)
        testSetPerturbations = [p for p in testSetEligible if p not in excessTestEligible]                      
        trainingSetPerturbations = list(testSetIneligible) + list(excessTestEligible) 
    # Now that the random part is done, we can start using sets. Order may change but content won't. 
    testSetPerturbations     = set(testSetPerturbations)
    trainingSetPerturbations = set(trainingSetPerturbations)
    adata_train    = adata[adata.obs["perturbation"].isin(trainingSetPerturbations),:]
    adata_heldout  = adata[adata.obs["perturbation"].isin(testSetPerturbations),    :]
    adata_train.uns[  "perturbed_and_measured_genes"]     = set(adata_train.uns[  "perturbed_and_measured_genes"]).intersection(trainingSetPerturbations)
    adata_heldout.uns["perturbed_and_measured_genes"]     = set(adata_heldout.uns["perturbed_and_measured_genes"]).intersection(testSetPerturbations)
    adata_train.uns[  "perturbed_but_not_measured_genes"] = set(adata_train.uns[  "perturbed_but_not_measured_genes"]).intersection(trainingSetPerturbations)
    adata_heldout.uns["perturbed_but_not_measured_genes"] = set(adata_heldout.uns["perturbed_but_not_measured_genes"]).intersection(testSetPerturbations)
    print("Test set size:")
    print(len(testSetPerturbations))
    print("Training set size:")
    print(len(trainingSetPerturbations))
    return adata_train, adata_heldout

def averageWithinPerturbation(ad: anndata.AnnData, confounders = []):
    """Average the expression levels within each level of ad.obs["perturbation"].

    Args:
        ad (anndata.AnnData): Object conforming to the validity checks in the load_perturbations module.
    """
    if len(confounders) != 0:
        raise NotImplementedError("Haven't yet decided how to handle confounders when merging replicates.")

    perts = ad.obs["perturbation"].unique()
    new_ad = anndata.AnnData(
        X = np.zeros((len(perts), len(ad.var_names))),
        obs = pd.DataFrame(
            {"perturbation":perts}, 
            index = perts, 
            columns=ad.obs.columns.copy(),
        ),
        var = ad.var,
    )
    for p in perts:
        p_idx = ad.obs["perturbation"]==p
        new_ad[p,].X = ad[p_idx,:].X.mean(0)
        new_ad.obs.loc[p,:] = ad[p_idx,:].obs.iloc[0,:]
        new_ad.obs.loc[p,"expression_level_after_perturbation"] = ad.obs.loc[p_idx, "expression_level_after_perturbation"].mean()
    new_ad.obs = new_ad.obs.astype(dtype = {c:ad.obs.dtypes[c] for c in new_ad.obs.columns}, copy = True)
    new_ad.raw = ad.copy()
    new_ad.uns = ad.uns.copy()
    return new_ad


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
