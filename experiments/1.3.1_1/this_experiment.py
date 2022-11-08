import seaborn as sns
import pandas as pd
import scanpy as sc
import os
import sys
import predict
import evaluator
import gc
import scanpy as sc

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
  # Redo clustering at different resolutions
  clusterResolutions = []
  num_clusters = []
  for r in [i/10.0 for i in range(1,11,2)]:
    sc.tl.leiden(train_data, resolution=r)
    new_name = "leiden_r="+str(round(r, 1))
    clusterResolutions.append(new_name)
    num_clusters.append(train_data.obs["leiden"].nunique())
    train_data.obs[new_name] = train_data.obs["leiden"]
    sc.pl.umap(train_data, color = "leiden")
    train_data.obs[new_name] = train_data.obs[new_name].astype("category")
    train_data.uns[new_name + "_colors"] = train_data.uns["leiden_colors"]
  n_resolutions = len(clusterResolutions)
  # Train on coarser or finer clusters

  experiments = pd.DataFrame(
    {
      "network":      ([n             for n in networks.keys()] + [list(networks.keys())[0]])*len(clusterResolutions),
      "network_prior":(["restrictive" for _ in networks.keys()] + ["ignore"]                )*len(clusterResolutions),
      "cluster_resolution":[f for f in clusterResolutions for _ in range(len(networks.keys()) + 1) ],
      "num_clusters":      [f for f in num_clusters       for _ in range(len(networks.keys()) + 1) ]
    }
  )
  predictions = {}
  for i in experiments.index:
    grn = predict.GRN(train=train_data, network=networks[experiments.loc[i,'network']])
    grn.extract_features(method = "tf_rna")
    grn.fit(
        method = "linear", 
        cell_type_labels = experiments.loc[i, "cluster_resolution"],
        cell_type_sharing_strategy = "distinct",
        network_prior = experiments.loc[i,'network_prior'],
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