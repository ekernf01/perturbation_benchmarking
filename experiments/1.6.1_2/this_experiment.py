import sys
import os
import numpy as np
import pandas as pd
import evaluator
import ggrn
import scanpy as sc
import anndata
import gc 
import torch

os.environ['WANDB_MODE'] = 'offline'
os.environ["WANDB_SILENT"] = "true"
import warnings
warnings.filterwarnings('ignore') # setting ignore as a parameter



def lay_out_runs(
  train_data: anndata.AnnData, 
  test_data: anndata.AnnData, 
  networks: dict, 
  outputs: str,
  metadata: dict,
) -> pd.DataFrame:
    """Lay out the specific runs done in this experiment.
    Args:
        train_data (anndata.AnnData): _description_
        test_data (anndata.AnnData):  usually not used, except in weird cases like the "oracle structure" experiment
        perturbationsToPredict (list):  genes and the expression level to set them to, e.g. {("FOXN1", 0), ("PAX9", 0)}
        networks (dict): dict with string keys and LightNetwork values
        outputs (str): folder name to save results in
        metadata (dict): metadata for this Experiment, from metadata.json. See this repo's global README.
    Returns:
        pd.DataFrame: metadata on the different conditions in this experiment
    """
    experiments = pd.DataFrame(
        {
            "regression_method":[
                "DCDFG-spectral_radius-linearlr-False",   # NOTEARS-LR
                "DCDFG-spectral_radius-mlplr-False",      # DCDFG
                "mean", 
                "median", 
                # "DCDFG-exp-linear-False",                 # NOTEARS
                # "DCDFG-exp-linear-True"                   # NOBEARS
            ]
        }
    )
    return experiments


def do_one_run(
  experiments: pd.DataFrame, 
  i: int, 
  train_data: anndata.AnnData, 
  test_data: anndata.AnnData, 
  networks: dict, 
  outputs: str,
  metadata: dict, 
  ) -> anndata.AnnData:
    """Do one run (fit a GRN model and make predictions) as part of this experiment.
    Args:
        experiments (pd.DataFrame): Output of lay_out_runs
        i (int): A value from the experiments.index
        Other args: see help(lay_out_runs)
    Returns:
        anndata.AnnData: Predicted expression
    """
    grn = ggrn.GRN(
        train=train_data,
    )
    grn.extract_tf_activity(method = "tf_rna")
    assert metadata["regression_method"] is None # We want to get this from the csv, not the metadata
    grn.fit(
        method = experiments.loc[i, "regression_method"], 
        cell_type_labels = None,
        cell_type_sharing_strategy = "identical",
        network_prior = "ignore",
        pruning_strategy = "none", 
        projection = "none", 
        kwargs = { 
            "num_train_epochs": 20000, 
            "num_fine_epochs" : 10000,
            "num_gpus": 1 if torch.cuda.is_available() else 0,
            "train_batch_size": 64
        }
    )
    return grn


def plot(evaluationResults, output):
    """Plots specific to this experiment.
    Args:
    evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
    output (_type_): where to save output
    """
