import gc
import seaborn as sns
import pandas as pd
import os
import sys
import scanpy as sc
# PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
# sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'benchmarking', 'src'))) 
import evaluator

def run(train_data, test_data, perturbationsToPredict, networks, outputs):
  """Prediction code specific to this experiment.
  Args:
    - train_data
    - test_data: usually not used, except in weird cases like the "oracle structure" experiment
    - perturbationsToPredict: genes and the expression level to set them to, e.g. {("FOXN1", 0), ("PAX9", 0)}
    - networks
    - outputs: folder name to save results in
  Return: 
    - experiments: metadata on the different conditions in this experiment
    - predictions: dictionary with keys equal to index of experiments. 
        Each value is an AnnData object, and its .obs must have the same index, "perturbation", and "expression_level_after_perturbation" as test_data.
    - other: can be anything
  """
  # One more network used in a reprogramming-related project is a union of FANTOM4 and STRING.
  networks["mogrify"] = pd.concat([networks[n] for n in ['MARA_FANTOM4','STRING']])
  n_networks = len(networks.keys())
  experiments = pd.DataFrame(
    {
      "network":[n for n in networks.keys()], 
    })
  predictions = {}
  for i in experiments.index:
    grn = predict.GRN(
      train=train_data, 
      network=networks[experiments.loc[i,'network']]
    )
    grn.extract_features(method = "tf_rna")
    grn.fit(
        method = "linear", 
        cell_type_labels = None,
        cell_type_sharing_strategy = "identical",
        network_prior = "restrictive",
        pruning_strategy = "none", 
        pruning_parameter = None,
        projection = "none", 
    )
    predictions[i] = grn.predict(perturbationsToPredict) 
    del grn
    gc.collect()
  other = None
  return experiments, predictions, other




def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """
  baseNetworkComparisonFigure = sns.FacetGrid(evaluationResults[~evaluationResults['somePredictionRefused']], 
                                              col = 'pruning',
                                              sharey = False, 
                                              height=5, 
                                              aspect=1).set(title = "Performance")
  baseNetworkComparisonFigure.map(sns.violinplot, "spearman", "network", 
                                  palette=["r", "b", "k", "y", "g"]
                              ).add_legend()
  baseNetworkComparisonFigure.set(ylabel="Spearman correlation\nminus average over all methods")
  baseNetworkComparisonFigureCompact.savefig(os.path.join(output, "baseNetworkComparisonFigure.pdf"))
  baseNetworkComparisonFigureCompact.savefig(os.path.join(output, "baseNetworkComparisonFigure.svg"))
  summary = evaluationResults[~evaluationResults['somePredictionRefused']]
  summary = summary.groupby(["pruning", "network"]).mean()[["spearman"]].reset_index(["pruning", "network"])
  summary = summary.merge(network_sizes)
  summary.sort_values(['pruning', 'network'], inplace=True)
  summary.to_csv(os.path.join(outputs, "networksExperimentEvaluationSummary.csv"))
  baseNetworkComparisonFigureCompact = sns.scatterplot(data=summary[[p!="harsh" for p in summary["pruning"]]],
                  x='numEdges',
                  y='spearman', 
                  hue='network')
  baseNetworkComparisonFigureCompact.set_xscale("log")
  baseNetworkComparisonFigureCompact.set(title="Density vs performance")
  baseNetworkComparisonFigureCompact.legend(loc='center left', bbox_to_anchor=(1, 0.5))
  baseNetworkComparisonFigureCompact.savefig(os.path.join(output, "baseNetworkComparisonFigureCompact.pdf"))
  baseNetworkComparisonFigureCompact.savefig(os.path.join(output, "baseNetworkComparisonFigureCompact.svg"))
