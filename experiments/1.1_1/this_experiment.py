import seaborn as sns
import pandas as pd
import sys
import os
import numpy as np
import evaluator

def run(train_data, test_data, perturbationsToPredict, networks, outputs):
  """Prediction code specific to this experiment.
  Args:
    - train_data
    - test_data
    - networks
    - outputs: folder name to save results in
  Return: 
    - experiments: metadata on the different conditions in this experiment
    - predictions: dictionary with keys equal to index of experiments. 
        Each value is an AnnData object, and its .obs must have the same index, "perturbation", and "expression_level_after_perturbation" as test_data.
    - other: can be anything
  """
  network_sizes = pd.DataFrame({bn:evaluator.countMatrixEdges(networks[bn]) for bn in networks}, index = ["numEdges"])
  network_sizes = network_sizes.T.reset_index().rename({"index":"network"}, axis = 1)

  threshold_number = [int(f) for f in np.logspace(np.log10(20000), np.log10(network_sizes['numEdges'].max()), 10)]
  n_pruning = len(threshold_number)
  experiments = pd.DataFrame({"threshold_number":threshold_number,
                              "network":["dense"]*n_pruning, 
                              "p":[1]*n_pruning,
                              "pruning":["none"]*n_pruning})

  experiments["log10_n_edges"] = round(np.log10(experiments["threshold_number"]), 2)
  experiments.to_csv(os.path.join(outputs, "networkExperiments.csv"))
  predictions = {
      i: evaluator.trainCausalModelAndPredict(
        expression=train_data,
        baseNetwork=networks[experiments.loc[i,'network']],
        memoizationName=os.path.join(outputs, str(i) + ".celloracle.oracle"), 
        perturbations=perturbationsToPredict,
        clusterColumnName = "fake_cluster",
        pruningParameters = {
          "p":experiments.loc[i,'p'], 
          "threshold_number":experiments.loc[i,'threshold_number']
          }
        ) 
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