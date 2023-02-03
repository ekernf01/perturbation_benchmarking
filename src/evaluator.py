"""evaluator.py is a collection of functions for making and testing predictions about expression fold change after genetic perturbations.
It dedicates particular attention to interfacing with CellOracle, a thought-provoking and flexible perturbation prediction method.
"""
from joblib import Parallel, delayed, cpu_count
import numpy as np
import pandas as pd
import anndata
from scipy.stats import spearmanr as spearmanr
import os 
import sys
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
    metrics = ['spearman', 'mse', 'mae', 'mae_benefit']
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
        vlnplot[metric].save(f'{outputs}/{metric}.svg', method = "selenium")
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
            right_on="gene")

    # Proteoform diversity information is not yet used because it would be hard to summarize this into a numeric measure of complexity.
    # But this code may be useful if we aim to continue that work later on.
    proteoform_diversity = pd.read_csv("../accessory_data/uniprot-compressed_true_download_true_fields_accession_2Cid_2Cprotei-2023.02.02-15.27.12.44.tsv.gz", sep = "\t")
    proteoform_diversity.head()
    proteoform_diversity_summary = pd.DataFrame(
        {
            "is_glycosylated": ~proteoform_diversity["Glycosylation"].isnull(),
            "has_ptm": ~proteoform_diversity["Post-translational modification"].isnull(),
        },
        index = proteoform_diversity.index,
    )
    proteoform_diversity_characteristics = proteoform_diversity_summary.columns.copy()

    # measures of evolutionary constraint 
    evolutionary_characteristics = ["pLI"]
    evolutionary_constraint = pd.read_csv("../accessory_data/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz", sep = "\t")
    evolutionary_constraint = evolutionary_constraint.groupby("gene").agg(func = max)
    if any(not x in df.columns for x in evolutionary_characteristics):
        df = pd.merge(
            evolutionary_constraint,
            df.copy(),
            left_on="gene", 
            right_on="gene")
    
    # measures of connectedness
    degree = pd.read_csv("../accessory_data/degree_info.csv.gz")
    degree = degree.rename({"Unnamed: 0":"gene"}, axis = 1)
    degree["gene"] = [str(g).upper() for g in degree["gene"]]
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
            right_on="gene")
    try:
        df.reset_index(inplace=True)
    except:
        pass
    return df, {
        "evolutionary_characteristics":evolutionary_characteristics,
        "expression_characteristics": expression_characteristics, 
        "degree_characteristics": degree_characteristics,
        }

def studyPredictableGenes(evaluationPerTarget, train_data, save_path, factor_varied, default_level, genes_considered_as):
    os.makedirs(os.path.join(save_path, genes_considered_as, "MAE_determinants"), exist_ok=True)
    # Plot various factors against our per-gene measure of predictability 
    evaluationPerTarget, types_of_gene_data = addGeneMetadata(evaluationPerTarget, train_data)
    types_of_gene_data["out-degree"] = [s for s in types_of_gene_data["degree_characteristics"] if "out-degree" in s]
    types_of_gene_data["in-degree"] = [s for s in types_of_gene_data["degree_characteristics"] if "in-degree" in s]
    for t in types_of_gene_data.keys():
        long_data = pd.melt(
            evaluationPerTarget, 
            id_vars=["model_beats_mean_on_this_gene", factor_varied], 
            value_vars=types_of_gene_data[t], 
            var_name='property_of_gene', 
            value_name='value', 
            col_level=None, 
            ignore_index=True)
        long_data[f"{factor_varied}_"] = long_data[factor_varied] + "__" + long_data["model_beats_mean_on_this_gene"].astype(str)
        try:
            long_data = long_data.query(f"{factor_varied} != '{default_level}'")
        except:
            pass
        chart = alt.Chart(long_data).mark_boxplot().encode(
                x = f'{factor_varied}_:N',
                y = "value:Q",
                color='model_beats_mean_on_this_gene:N',
            ).facet(
                "property_of_gene:N", 
                columns= 10,
            ).resolve_scale(
                y='independent'
            )
        _ = alt.data_transformers.disable_max_rows()
        chart.save(os.path.join(save_path, genes_considered_as, f"predictability_vs_{t}.svg"), method = "selenium")

    # How many genes are we just predicting a constant for?
    if genes_considered_as == "targets":
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
            os.makedirs(os.path.join(save_path, genes_considered_as, "variety_in_predictions"), exist_ok=True)
            chart.save( os.path.join(save_path, genes_considered_as, "variety_in_predictions", f"{condition}.svg"), method = "selenium")

    # Gene set enrichments on best-predicted genes
    for condition in evaluationPerTarget[factor_varied].unique():
        os.makedirs(os.path.join(save_path, genes_considered_as, "enrichr_on_best", condition), exist_ok=True)
        gl = evaluationPerTarget.loc[evaluationPerTarget[factor_varied]==condition]
        gl = list(gl.sort_values("mae_benefit", ascending=False).head(50)["gene"].unique())
        pd.DataFrame(gl).to_csv(os.path.join(save_path, genes_considered_as, "enrichr_on_best", condition, "input_genes.txt"), index = False,  header=False)
        for gene_sets in ['GO Molecular Function 2021', 'GO Biological Process 2021', 'Jensen TISSUES', 'ARCHS4 Tissues', 'Chromosome Location hg19']:
            try:
                _ = gseapy.enrichr(
                    gene_list=gl,
                    gene_sets=gene_sets.replace(" ", "_"), 
                    outdir=os.path.join(save_path, genes_considered_as, "enrichr_on_best", condition, f"{gene_sets}"), 
                    format='svg',
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
    ).save(os.path.join(outputs, gene + ".svg"), method = "selenium")
    return   

def postprocessEvaluations(evaluations, experiments, factor_varied, default_level):
    evaluations   = pd.concat(evaluations)
    evaluations   = evaluations.merge(experiments,   how = "left", right_index = True, left_on = "index")
    evaluations   = pd.DataFrame(evaluations.to_dict())
    # Add some info on each evaluation-per-target, such as the baseline MAE
    evaluations["target"] = [i[1] for i in evaluations.index]
    is_baseline = [f==default_level for f in evaluations[factor_varied]]
    evaluations["mae_baseline"] = np.NaN
    evaluations.loc[is_baseline, "mae_baseline"] = evaluations.loc[is_baseline, "mae"]
    evaluations = evaluations.groupby("target", group_keys=False).apply(
        lambda x:
            x.fillna(
                x.loc[x[factor_varied] == default_level, "mae"].values[0]
            )
    )
    evaluations["mae_benefit"] = evaluations["mae_baseline"] - evaluations["mae"]
    evaluations = evaluations.sort_values("mae_benefit", ascending=False)
    evaluations["model_beats_mean_on_this_gene"] = evaluations["mae_benefit"]>0
    # Sometimes these are processed by the same code downstream and it's convenient to have a "gene" column.
    try:
        evaluations["gene"] = evaluations["target"]
    except KeyError:
        pass
    try:
        evaluations["gene"] = evaluations["perturbation"]
    except KeyError:
        pass
    return evaluations

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
    if default_level is None:
        default_level = experiments.loc[0,factor_varied]
    evaluationPerPert = {}
    evaluationPerTarget = {}
    shared_var_names = list(set.intersection(*[set(predictions[experiment].var_names) for experiment in predictions.keys()]))    
    evaluations = Parallel(n_jobs=cpu_count()-1, verbose = 1, backend="loky")(
        delayed(evaluateOnePrediction)(
            expression =            heldout[experiment][:, shared_var_names],
            predictedExpression=predictions[experiment][:, shared_var_names],
            baseline =             baseline[experiment][:, shared_var_names],
            doPlots=do_scatterplots,
            outputs = outputs,
            experiment_name = experiment,
            classifier=classifier,        
        )
        for experiment in predictions.keys()
    )
    # That parallel code returns a list of tuples. I want a pair of dicts instead. 
    for i,experiment in enumerate(predictions.keys()):
        evaluationPerPert[experiment], evaluationPerTarget[experiment] = evaluations[i]
        evaluationPerPert[experiment]["index"]   = experiment
        evaluationPerTarget[experiment]["index"] = experiment
    del evaluations
    # Concatenate and add some extra info
    evaluationPerPert = postprocessEvaluations(evaluationPerPert, experiments, factor_varied, default_level)
    evaluationPerTarget = postprocessEvaluations(evaluationPerTarget, experiments, factor_varied, default_level)
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
    perturbation_plot_path = os.path.join(outputs, "perturbations", str(experiment_name))
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
            scatterplot.save(os.path.join(perturbation_plot_path, f"{pert}.svg"), method = "selenium")
            if is_easiest:
                scatterplot.save(os.path.join(perturbation_plot_path, f"_easiest({pert}).svg"), method = "selenium")
            if is_hardest:
                scatterplot.save(os.path.join(perturbation_plot_path, f"_hardest({pert}).svg"), method = "selenium")
    metrics["perturbation"] = metrics.index
    return metrics, metrics_per_target
    