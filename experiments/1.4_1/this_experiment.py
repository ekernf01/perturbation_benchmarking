import seaborn as sns
import pandas as pd
import os
import sys
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
      "p":[1]*n_networks,
      "threshold_number":[int(network_sizes['numEdges'].max())]*n_networks,
      "pruning":["none"]*n_networks
    })
  predictions = {
    i: evaluator.trainCausalModelAndPredict(
      expression=train_data,
      baseNetwork=networks[experiments.loc[i,'network']],
      memoizationName=os.path.join(outputs, str(i), ".celloracle.oracle"),
      perturbations=perturbationsToPredict,
      clusterColumnName = "fake_cluster",
      pruningParameters = {
        "p":experiments.loc[i,'p'], 
        "threshold_number":experiments.loc[i,'threshold_number']
      }) 
    for i in experiments.index
  }
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
