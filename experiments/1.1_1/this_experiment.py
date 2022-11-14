import sys
import os
import numpy as np
import pandas as pd
import evaluator
import predict
import scanpy as sc
import anndata
import gc 

def lay_out_runs(
  train_data: anndata.AnnData, 
  test_data: anndata.AnnData, 
  perturbationsToPredict: list, 
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
  size_of_dense_network = len(train_data.var_names)*len(grn.tf_list)
  threshold_number = [int(f) for f in np.logspace(np.log10(20000), np.log10(size_of_dense_network), 10)]
  experiments = pd.DataFrame({"threshold_number":threshold_number})
  experiments["log10_n_edges"] = round(np.log10(experiments["threshold_number"]), 2)
  return experiments
  
def do_one_run(
  experiments: pd.DataFrame, 
  i: int, 
  train_data: anndata.AnnData, 
  test_data: anndata.AnnData, 
  perturbationsToPredict: list, 
  networks: dict, 
  outputs: str
) -> anndata.AnnData:
  print("Running setting " + str(i))
  sys.stdout.flush()
  grn = predict.GRN(train=train_data)
  grn.extract_features(method = "tf_rna")
  grn.fit(
      method = "linear", 
      cell_type_labels = None,
      cell_type_sharing_strategy = "identical",
      network_prior = "ignore",
      pruning_strategy = "prune_and_refit", 
      pruning_parameter = experiments.loc[i,'threshold_number'],
      projection = "none", 
  )
  predictions = grn.predict(perturbationsToPredict) 
  del grn
  gc.collect()
  return predictions

def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """