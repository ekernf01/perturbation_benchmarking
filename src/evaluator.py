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
import gseapy
sys.path.append("src")
import altair as alt

def makeMainPlots(
    evaluationPerPert: pd.DataFrame, 
    evaluationPerTarget: pd.DataFrame, 
    outputs: str, 
    factor_varied:str, 
    facet_by: str = None, 
    color_by: str = None, 
    metrics = ['spearman', 'mse', 'mae']
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
        if color_by is not None:
            evaluationPerPert[factor_varied] = [str(a) + str(b) for a,b in zip(evaluationPerPert[factor_varied], evaluationPerPert[color_by])]
        if facet_by is not None:
            group_mean_by.append(facet_by)
        means = evaluationPerPert.groupby(group_mean_by, as_index=False)[[metric]].mean()
        vlnplot[metric] = alt.Chart(
                data = evaluationPerPert, 
                title = f"{metric} (predicted log fold change vs observed)"
            ).mark_boxplot(extent='min-max')
        # Faceting fights with layering, so skip the means if faceting.
        if facet_by is None:
            vlnplot[metric] = vlnplot[metric] + alt.Chart(data = means).mark_point(color="black")
        if color_by is not None:
            vlnplot[metric]=vlnplot[metric].encode(
                y=alt.Y(f'{metric}:Q'),
                color=color_by + ':N',
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
                columns=int(np.ceil(np.sqrt(len(evaluationPerPert[facet_by].unique())))), 
            )
        vlnplot[metric].save(f'{outputs}/{metric}.html')
    return vlnplot

def addGeneMetadata(df, adata):

    # Measures derived from the expression data, e.g. overdispersion
    expression_characteristics = [
        'highly_variable', 'highly_variable_rank', 'means',
        'variances', 'variances_norm'
    ]
    if any(not x in df.columns for x in expression_characteristics):
        df = pd.merge(
            adata.var,
            df.copy(),
            left_index=True, 
            right_on="target")

    # measures of evolutionary constraint 
    evolutionary_characteristics = ["pLI"]
    evolutionary_constraint = pd.read_csv("../accessory_data/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz", sep = "\t")
    evolutionary_constraint = evolutionary_constraint.groupby("gene").agg(func = max)
    if any(not x in df.columns for x in evolutionary_characteristics):
        df = pd.merge(
            evolutionary_constraint,
            df.copy(),
            left_on="gene", 
            right_on="target")
    
    # measures of connectedness
    degree = pd.read_csv("../accessory_data/degree_info.csv.gz")
    degree = degree.rename({"Unnamed: 0":"gene"}, axis = 1)
    degree = degree.pivot_table(
        index=['gene'], 
        values=['in-degree', 'out-degree'], 
        columns=['network']
    )
    degree.fillna(0)
    degree.columns = ['_'.join(col) for col in degree.columns.values]
    degree_characteristics = list(degree.columns)
    if any(not x in df.columns for x in degree_characteristics):
        df = pd.merge(
            degree,
            df.copy(),
            left_on="gene", 
            right_on="target")
    try:
        df.reset_index(inplace=True)
    except:
        pass
    types_of_gene_data = evolutionary_characteristics + expression_characteristics + degree_characteristics
    return df, types_of_gene_data

def studyPredictableGenes(evaluationPerTarget, train_data, save_path, factor_varied):
    os.makedirs(os.path.join(save_path, "target_genes_best_worst", "MAE_determinants"), exist_ok=True)
    # Plot various factors against our per-gene measure of predictability 
    evaluationPerTarget, types_of_gene_data = addGeneMetadata(evaluationPerTarget, train_data)
    for x in types_of_gene_data:
        chart = alt.Chart(evaluationPerTarget).mark_bar().encode(
            x=alt.X(f"{x}:Q", bin=True),
            y=alt.Y('count()', stack='normalize'),
            color='model_beats_mean_on_this_gene',
        ).facet(
            factor_varied, 
            columns = 4,
        ).properties(
            title=f'{x} versus modeling success'
        )
        _ = alt.data_transformers.disable_max_rows()
        chart.save(os.path.join(save_path, "target_genes_best_worst", f"MAE_determinants/{x}_.html"))
        chart = alt.Chart(evaluationPerTarget).mark_bar().encode(
            x=alt.X(f"{x}:Q", bin=True),
            y=alt.Y('count()'),
            color='model_beats_mean_on_this_gene',
        ).facet(
            factor_varied, 
            columns = 4,
        ).properties(
            title=f'{x} versus modeling success'
        )
        _ = alt.data_transformers.disable_max_rows()
        chart.save(os.path.join(save_path, "target_genes_best_worst", f"MAE_determinants/{x}.html"))
    
    # Gene set enrichments on best-predicted targets
    cutoff = 0.01
    for condition in evaluationPerTarget[factor_varied].unique():
        subset = evaluationPerTarget.loc[evaluationPerTarget[factor_varied]==condition]
        n_constant = (subset["standard_deviation"] < cutoff).sum()
        n_total = subset.shape[0]  
        chart = alt.Chart(subset).mark_bar().encode(
                x=alt.X("standard_deviation:Q", bin=alt.BinParams(maxbins=30), scale=alt.Scale(type="sqrt")),
                y=alt.Y('count()'),
            ).properties(
                title=f'Standard deviation of predictions ({n_constant}/{n_total} are within {cutoff} of 0)'
            )
        _ = alt.data_transformers.disable_max_rows()
        os.makedirs(os.path.join(save_path, "target_genes_best_worst", "variety_in_predictions"), exist_ok=True)
        chart.save( os.path.join(save_path, "target_genes_best_worst", "variety_in_predictions", f"{condition}.html"))

    for condition in evaluationPerTarget[factor_varied].unique():
        os.makedirs(os.path.join(save_path, "target_genes_best_worst", "enrichr_on_best", condition), exist_ok=True)
        gl = evaluationPerTarget.loc[evaluationPerTarget[factor_varied]==condition]
        gl = list(gl.sort_values("mae_benefit", ascending=False).head(50)["Symbol"].unique())
        pd.DataFrame(gl).to_csv(os.path.join(save_path, "target_genes_best_worst", "enrichr_on_best", condition, "input_genes.txt"), index = False,  header=False)
        for gene_sets in ['GO Molecular Function 2021', 'GO Biological Process 2021', 'Jensen TISSUES', 'ARCHS4 Tissues', 'Chromosome Location hg19']:
            try:
                _ = gseapy.enrichr(
                    gene_list=gl,
                    gene_sets=gene_sets.replace(" ", "_"), 
                    outdir=os.path.join(save_path, "target_genes_best_worst", "enrichr_on_best", condition, f"{gene_sets}"), 
                    format='png',
                )
            except ValueError:
                pass
    return evaluationPerTarget

def plotOneTargetGene(gene, outputs, experiments, factor_varied, train_data, heldout_data, fitted_values, predictions):
    """For one gene, plot predicted + observed values for train + test."""
    expression = {
        e:pd.DataFrame({
            "index": [i for i in range(
                fitted_values[e][:,gene].shape[0] + 
                predictions[e][:,gene].shape[0]
            )],
            "experiment": e,
            "observed": np.concatenate([
                train_data[e][:,gene].X.squeeze(), 
                heldout_data[e][:,gene].X.squeeze(), 
            ]), 
            "predicted": np.concatenate([
                fitted_values[e][:,gene].X.squeeze(), 
                predictions[e][:,gene].X.squeeze(), 
            ]), 
            "is_trainset": np.concatenate([
                np.ones (fitted_values[e][:,gene].shape[0]), 
                np.zeros(  predictions[e][:,gene].shape[0]), 
            ]), 
        }) for e in predictions.keys() 
    }
    expression = pd.concat(expression)
    expression = expression.reset_index()
    expression = expression.merge(experiments, left_on="experiment", right_index=True)
    os.makedirs(os.path.join(outputs), exist_ok=True)
    alt.Chart(data=expression).mark_point().encode(
        x = "observed:Q",y = "predicted:Q", color = "is_trainset:N"
    ).properties(
        title=gene
    ).facet(
        facet = factor_varied, 
        columns=3,
    ).save(os.path.join(outputs, gene + ".html"))
    return   

def evaluateCausalModel(
    heldout:dict, 
    predictions:dict, 
    baseline:dict,
    experiments: pd.DataFrame, 
    outputs: str, 
    factor_varied: str,
    default_level = None, 
    classifier = None, 
    do_scatterplots = True):
    """Compile plots and tables comparing heldout data and predictions for same. 

    Args:
        heldout, predictions, baseline: each of these is a dictionary with keys equal to index of experiments. 
            Each value is an AnnData object. 
            Baseline is expression before perturbation, for use in calculating log fold change. 
        classifier (sklearn.LogisticRegression): Optional, to judge results on cell type accuracy. 
        experiments (pd.DataFrame): Metadata for the different combinations used in this experiment. 
        default_level: a value found in experiments.loc[:,factor_varied] to compare all others against. This is for use in repeated-measures analyses to remove the variation from gene to gene. 
        outputs (String): Saves output here.
    """
    # Get spearman and classifier accuracy 
    if default_level is None:
        default_level = experiments.loc[0,factor_varied]
    evaluationPerPert = {}
    evaluationPerTarget = {}
    shared_var_names = list(set.intersection(*[set(predictions[experiment].var_names) for experiment in predictions.keys()]))
    for experiment in predictions.keys(): 
        evaluationPerPert[experiment], evaluationPerTarget[experiment] = \
            evaluateOnePrediction(
                expression =            heldout[experiment][:, shared_var_names],
                predictedExpression=predictions[experiment][:, shared_var_names],
                baseline =             baseline[experiment][:, shared_var_names],
                doPlots=do_scatterplots,
                outputs = outputs,
                experiment_name = experiment,
                classifier=classifier,
            )
        evaluationPerPert[experiment]["index"]   = experiment
        evaluationPerTarget[experiment]["index"] = experiment
    evaluationPerPert   = pd.concat(evaluationPerPert)
    evaluationPerTarget = pd.concat(evaluationPerTarget)
    evaluationPerPert   = evaluationPerPert.merge(experiments,   how = "left", right_index = True, left_on = "index")
    evaluationPerTarget = evaluationPerTarget.merge(experiments, how = "left", right_index = True, left_on = "index")
    evaluationPerPert   = pd.DataFrame(evaluationPerPert.to_dict())
    evaluationPerTarget = pd.DataFrame(evaluationPerTarget.to_dict())
    # Add some info on each evaluation-per-target, such as the baseline MAE
    evaluationPerTarget["target"] = [i[1] for i in evaluationPerTarget.index]
    is_baseline = [f==default_level for f in evaluationPerTarget[factor_varied]]
    evaluationPerTarget["mae_baseline"] = np.NaN
    evaluationPerTarget.loc[is_baseline, "mae_baseline"] = evaluationPerTarget.loc[is_baseline, "mae"]
    evaluationPerTarget = evaluationPerTarget.groupby("target", group_keys=False).apply(
        lambda x:
            x.fillna(
                x.loc[x[factor_varied] == default_level, "mae"].values[0]
            )
    )
    evaluationPerTarget["mae_benefit"] = evaluationPerTarget["mae_baseline"] - evaluationPerTarget["mae"]
    evaluationPerTarget = evaluationPerTarget.sort_values("mae_benefit", ascending=False)
    evaluationPerTarget["model_beats_mean_on_this_gene"] = evaluationPerTarget["mae_benefit"]>0
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
    experiment_name: str,
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
                    outputs (str): Folder to save output in
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
    if not all(predictedExpression.obs_names == expression.obs_names):
        raise ValueError("expression and predictedExpression must have the same indices.")
    baseline = baseline.X.mean(axis=0).squeeze()
    metrics = pd.DataFrame(index = predictedExpression.obs.index, columns = ["spearman", "spearmanp", "cell_fate_correct", "mse"])
    metrics_per_target = pd.DataFrame(index = predictedExpression.var.index, columns = ["mae", "mse", "standard_deviation"])
    for target in predictedExpression.var.index:
        observed  = expression[         :,target].X.squeeze()
        predicted = predictedExpression[:,target].X.squeeze()
        metrics_per_target.loc[target,["standard_deviation"]] = np.std(predicted)
        metrics_per_target.loc[target,["mae"]] = np.abs(observed - predicted).sum()
        metrics_per_target.loc[target,["mse"]] = np.linalg.norm(observed - predicted)**2
    for pert in predictedExpression.obs.index:
        if do_careful_checks:
            if not all(
                                 expression.obs.loc[pert, ["perturbation", "expression_level_after_perturbation"]].fillna(0) == \
                        predictedExpression.obs.loc[pert, ["perturbation", "expression_level_after_perturbation"]].fillna(0) 
                    ):
                print(expression.obs.head()[["perturbation", "expression_level_after_perturbation"]])
                print(predictedExpression.obs.head()[["perturbation", "expression_level_after_perturbation"]])
                raise ValueError(f"Expression and predicted expression are different sizes or are differently name in experiment {experiment_name}.")
        observed  = expression[         pert,:].X.squeeze()
        predicted = predictedExpression[pert,:].X.squeeze()
        def is_constant(x):
            return np.std(x)<1e-12
        if type(predicted) is float and np.isnan(predicted) or is_constant(predicted - baseline) or is_constant(observed - baseline):
            metrics.loc[pert,["spearman","spearmanp", "cell_fate_correct", "mse", "mae", "proportion_correct_direction"]] = 0,1,np.nan,np.nan,np.nan,np.nan
        else:
            metrics.loc[pert,["spearman","spearmanp"]] = [x for x in spearmanr(observed - baseline, predicted - baseline)]
            metrics.loc[pert,"mse"] = np.linalg.norm(observed - predicted)**2
            metrics.loc[pert,"mae"] = np.abs(observed - predicted).mean()
            metrics.loc[pert,"proportion_correct_direction"] = np.mean((observed>=0) == (predicted >= 0))
            if classifier is not None:
                class_observed  = classifier.predict(np.reshape(observed,  (1, -1)))[0]
                class_predicted = classifier.predict(np.reshape(predicted, (1, -1)))[0]
                metrics.loc[pert,"cell_fate_correct"] = 1.0*(class_observed==class_predicted)  
    
    metrics["spearman"] = metrics["spearman"].astype(float)
    hardest = metrics["spearman"].idxmin()
    easiest = metrics["spearman"].idxmax()
    perturbation_plot_path = os.path.join(outputs, "samples_best_worst", str(experiment_name))
    for pert in metrics.index:
        observed  = expression[         pert,:].X.squeeze()
        predicted = predictedExpression[pert,:].X.squeeze()
        is_hardest = hardest==pert
        is_easiest = easiest==pert
        if doPlots | is_hardest | is_easiest:
            os.makedirs(perturbation_plot_path, exist_ok = True)
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
            alt.data_transformers.disable_max_rows()
            pd.DataFrame().to_csv(os.path.join(perturbation_plot_path, f"{pert}.txt"))
            scatterplot.save(os.path.join(perturbation_plot_path, f"{pert}.html"))
            if is_easiest:
                scatterplot.save(os.path.join(perturbation_plot_path, f"_easiest({pert}).html"))
            if is_hardest:
                scatterplot.save(os.path.join(perturbation_plot_path, f"_hardest({pert}).html"))
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


def splitData(adata, allowedRegulators, desired_heldout_fraction, type_of_split, data_split_seed):
    """Determine a train-test split satisfying constraints imposed by base networks and available data.
    
    A few factors complicate the training-test split. 

    - Perturbed genes may be absent from most base GRN's due to lack of motif information or ChIP data. 
        These perhaps should be excluded from the test data to avoid obvious failure cases.
    - Perturbed genes may not be measured. These perhaps should be excluded from test data because we can't
        reasonably separate their direct vs indirect effects.

    If type_of_split=="simple", we make no provision for dealing with the above concerns. The only restriction is that
    all controls go in the training set.
    If type_of_split=="interventional", the `allowedRegulators` arg can be specified in order to keep any user-specified
    problem cases out of the test data. No matter what, we still use those perturbed profiles as training data, hoping 
    they will provide useful info about attainable cell states and downstream causal effects. 

    For some collections of base networks, there are many factors ineligible for use as test data -- so many that 
    we use all the eligible ones for test and the only ineligible ones for training. 
    For other cases, such as dense base networks, we have more flexibility, so we send some perturbations to the 
    training set at random even if we would be able to use them in the test set.

    parameters:

    - adata (anndata.AnnData): Object satisfying the expectations outlined in the accompanying collection of perturbation data.
    - allowedRegulators (list or set): interventions allowed to be in the test set. 
    - type_of_split (str): if "interventional" (default), then any perturbation is placed in either the training or the test set, but not both. 
        If "simple", then we use a simple random split, and replicates of the same perturbation are allowed to go into different folds.

    """
    if desired_heldout_fraction is None:
        desired_heldout_fraction = 0.5
    if data_split_seed is None:
        data_split_seed = 0
    # For a deterministic result when downsampling an iterable, setting a seed alone is not enough.
    # Must also avoid the use of sets. 
    if type_of_split is None or np.isnan(type_of_split) or type_of_split == "interventional":
        get_unique_keep_order = lambda x: list(dict.fromkeys(x))
        allowedRegulators = [p for p in allowedRegulators if p in adata.uns["perturbed_and_measured_genes"]]
        testSetEligible   = [p for p in adata.obs["perturbation"] if     all(g in allowedRegulators for g in p.split(","))]
        testSetIneligible = [p for p in adata.obs["perturbation"] if not all(g in allowedRegulators for g in p.split(","))]
        allowedRegulators = get_unique_keep_order(allowedRegulators)
        testSetEligible   = get_unique_keep_order(testSetEligible)
        testSetIneligible = get_unique_keep_order(testSetIneligible)
        total_num_perts = len(testSetEligible) + len(testSetIneligible)
        eligible_heldout_fraction = len(testSetEligible)/(0.0+total_num_perts)
        if eligible_heldout_fraction < desired_heldout_fraction:
            print("Not enough profiles for the desired_heldout_fraction. Will use all available.")
            testSetPerturbations = testSetEligible
            trainingSetPerturbations = testSetIneligible
        elif eligible_heldout_fraction == desired_heldout_fraction: #nailed it
            testSetPerturbations = testSetEligible
            trainingSetPerturbations = testSetIneligible
        else:
            # Plenty of perts work for either.
            # Put some back in trainset to get the right size, even though we could use them in test set.
            numExcessTestEligible = int(np.ceil((eligible_heldout_fraction - desired_heldout_fraction)*total_num_perts))
            excessTestEligible = np.random.default_rng(seed=data_split_seed).choice(
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
    elif type_of_split == "simple":
        np.random.seed(data_split_seed)
        train_obs = np.random.choice(
            replace=False, 
            a = adata.obs_names, 
            size = round(adata.shape[0]*(1-desired_heldout_fraction)), 
        )
        for o in adata.obs_names:
            if adata.obs.loc[o, "is_control"]:
                train_obs = np.append(train_obs, o)
        test_obs = [i for i in adata.obs_names if i not in train_obs]
        adata_train    = adata[train_obs,:]
        adata_heldout  = adata[test_obs,:]
        trainingSetPerturbations = set(  adata_train.obs["perturbation"].unique())
        testSetPerturbations     = set(adata_heldout.obs["perturbation"].unique())
        adata_train.uns[  "perturbed_and_measured_genes"]     = set(adata_train.uns[  "perturbed_and_measured_genes"]).intersection(trainingSetPerturbations)
        adata_heldout.uns["perturbed_and_measured_genes"]     = set(adata_heldout.uns["perturbed_and_measured_genes"]).intersection(testSetPerturbations)
        adata_train.uns[  "perturbed_but_not_measured_genes"] = set(adata_train.uns[  "perturbed_but_not_measured_genes"]).intersection(trainingSetPerturbations)
        adata_heldout.uns["perturbed_but_not_measured_genes"] = set(adata_heldout.uns["perturbed_but_not_measured_genes"]).intersection(testSetPerturbations)
    else:
        raise ValueError(f"`type_of_split` must be 'simple' or 'interventional'; got {type_of_split}.")
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
        dtype = np.float32
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
