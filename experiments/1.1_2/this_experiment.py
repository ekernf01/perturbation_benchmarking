from turtle import down
import seaborn as sns
import pandas as pd
import sys
import os
import gc
import numpy as np
import scanpy as sc
import anndata
#sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'benchmarking', 'src'))) 
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
  downSampleFactors = [i/10.0 for i in range(3,11,2)]
  experiments = pd.DataFrame(
    {
      "network":      ([n             for n in networks.keys()] + [list(networks.keys())[0]])*len(downSampleFactors),
      "network_prior":(["restrictive" for _ in networks.keys()] + ["ignore"]                )*len(downSampleFactors),
      "training_set_size":[f for f in downSampleFactors for _ in range(len(networks.keys()) + 1) ]
    }
  )
  experiments = experiments.merge(pd.DataFrame({"seed":range(3)}), how='cross')
  predictions = {}
  for i in experiments.index:
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
    predictions[i] = grn.predict(perturbationsToPredict)  
    del grn
    gc.collect()
  other = None
  return experiments, predictions, other



def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """