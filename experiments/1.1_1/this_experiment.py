import seaborn as sns
import pandas as pd
import sys
import os
import numpy as np
import evaluator
import predict

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
  grn = predict.GRN(train=train_data)
  grn.extract_features(method = "tf_rna")
  print(len(train_data.var_names))
  print(len(grn.tf_list))
  size_of_dense_network = len(train_data.var_names)*len(grn.tf_list)
  threshold_number = [int(f) for f in np.logspace(np.log10(20000), np.log10(size_of_dense_network), 10)]
  experiments = pd.DataFrame({"threshold_number":threshold_number})
  experiments["log10_n_edges"] = round(np.log10(experiments["threshold_number"]), 2)
  experiments.to_csv(os.path.join(outputs, "networkExperiments.csv"))

  predictions = {}
  for i in experiments.index:
    grn.fit(
        method = "linear", 
        cell_type_labels = None,
        cell_type_sharing_strategy = "identical",
        network_prior = "ignore",
        pruning_strategy = "prune_and_refit", 
        pruning_parameter = experiments.loc[i,'threshold_number'],
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