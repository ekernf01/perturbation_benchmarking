import sys
import os
import numpy as np
import pandas as pd
import evaluator
import predict
import scanpy as sc
import anndata
import gc 
import altair as alt

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
  n_networks = len(networks.keys())
  experiments = pd.DataFrame(
    {
      "network":[n for n in networks.keys()], 
      "network_prior":[
        "ignore" if list(networks.keys())[i] == "dense" else "restrictive"
        for i in range(len(networks.keys()))
      ]
    }
  )
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
      network_prior = experiments.loc[i, 'network_prior'],
      pruning_strategy = "none", 
      pruning_parameter = None,
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
  vlnplot = alt.Chart(
      data = evaluationResults, 
      title = "Spearman correlation (predicted log fold change vs observed)"
  ).mark_boxplot()
  vlnplot=vlnplot.encode(
      y=alt.Y('spearman:Q'),
      color="is_tf_in_network:N",
      x=alt.X(
          'network:N'
      )
  ).properties(
      width=400,
      height=400
  )
  vlnplot.save(f'{output}/custom_stripchart.pdf')
  return vlnplot