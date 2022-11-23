import seaborn as sns
import pandas as pd
import sys
import os
import gc
import numpy as np
import evaluator
import predict
import getsize
import scanpy as sc
import psutil 
from memory_profiler import profile
import anndata

def lay_out_runs(
  train_data: anndata.AnnData, 
  test_data: anndata.AnnData, 
  networks: dict, 
  outputs: str
  ) -> pd.DataFrame:
  """Lay out the specific runs done in this experiment.

  Args:
      train_data (anndata.AnnData): _description_
      test_data (anndata.AnnData):  usually not used, except in weird cases like the "oracle structure" experiment
      perturbationsToPredict (list):  genes and the expression level to set them to, e.g. {("FOXN1", 0), ("PAX9", 0)}
      networks (dict): dict with string keys and LightNetwork values
      outputs (str): folder name to save results in

  Returns:
      pd.DataFrame: metadata on the different conditions in this experiment

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
  return experiments
  
def do_one_run(
  experiments: pd.DataFrame, 
  i: int, 
  train_data: anndata.AnnData, 
  test_data: anndata.AnnData, 
  networks: dict, 
  outputs: str
  ) -> anndata.AnnData:
  """Do one run (fit a GRN model and make predictions) as part of this experiment.

  Args:
      experiments (pd.DataFrame): Output of lay_out_runs
      i (int): A value from the experiments.index
      Other args: see help(lay_out_runs)

  Returns:
      anndata.AnnData: Predicted expression
  """
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
      do_parallel = True,
  )
  return grn


def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """