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
  downSampleFactors = [i/10.0 for i in range(3,11,2)]
  experiments = pd.DataFrame(
    {
      "network":      ([n             for n in networks.keys()] + [list(networks.keys())[0]])*len(downSampleFactors),
      "network_prior":(["restrictive" for _ in networks.keys()] + ["ignore"]                )*len(downSampleFactors),
      "training_set_size":[f for f in downSampleFactors for _ in range(len(networks.keys()) + 1) ]
    }
  )
  experiments = experiments.merge(pd.DataFrame({"seed":range(3)}), how='cross')
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
    train=evaluator.downsample(
        train_data,
        experiments.loc[i, "training_set_size"], 
        seed=experiments.loc[i, "seed"], 
      ), 
    network=networks[experiments.loc[i,'network']]
  )
  grn.extract_features(method = "tf_rna")
  grn.fit(
      method = "linear", 
      cell_type_labels = None,
      cell_type_sharing_strategy = "identical",
      network_prior = experiments.loc[i, "network_prior"],
      pruning_strategy = "none", 
      projection = "none", 
  )
  return grn



def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """