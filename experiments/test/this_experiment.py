import seaborn as sns
import pandas as pd
import sys
import os
import gc
import numpy as np
import evaluator
import predict
import scanpy as sc
import psutil 

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
  network_sizes = pd.DataFrame({bn:networks[bn].get_num_edges() for bn in networks}, index = ["numEdges"])
  network_sizes = network_sizes.T.reset_index().rename({"index":"network"}, axis = 1)

  threshold_number = [int(f) for f in np.logspace(np.log10(20000), np.log10(network_sizes['numEdges'].max()), 1)]
  n_pruning = len(threshold_number)
  experiments = pd.DataFrame({"threshold_number":threshold_number,
                              "network":["celloracle_human all"]*n_pruning, 
                              "p":[1]*n_pruning,
                              "pruning":["none"]*n_pruning})

  experiments["log10_n_edges"] = round(np.log10(experiments["threshold_number"]), 2)
  experiments.to_csv(os.path.join(outputs, "networkExperiments.csv"))
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
        pruning_strategy = "prune_and_refit", 
        pruning_parameter = experiments.loc[i,'threshold_number'],
        projection = "none", 
        do_parallel = False
    )
    predictions[i] = grn.predict(perturbationsToPredict)  
    print("Trying to deallocate memory.")
    print(psutil.Process().memory_info().rss)
    del grn
    print(gc.collect())
    print(psutil.Process().memory_info().rss)
  other = None
  return experiments, predictions, other



def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """