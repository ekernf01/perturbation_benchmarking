import seaborn as sns
import pandas as pd
import sys
import os
import evaluator
import predict

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
  n_networks = len(networks.keys())
  network_sizes = pd.DataFrame({bn:evaluator.countMatrixEdges(networks[bn]) for bn in networks}, index = ["numEdges"])
  network_sizes = network_sizes.T.reset_index().rename({"index":"network"}, axis = 1)

  experiments = pd.DataFrame({"network":[n for n in networks.keys()]})

  predictions = {}
  for i in experiments.index:
    grn = predict.GRN(train=train_data, network=networks[experiments.loc[i,'network']])
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

  other = None
  return experiments, predictions, other



def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """