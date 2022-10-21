import seaborn as sns
import pandas as pd
import scanpy as sc
import os
import sys
import evaluator

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
  n_networks = len(networks)
  network_sizes = pd.DataFrame({bn:evaluator.countMatrixEdges(networks[bn]) for bn in networks}, index = ["numEdges"])
  network_sizes = network_sizes.T.reset_index().rename({"index":"network"}, axis = 1)
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
  experiments = pd.DataFrame({"network":[n for n in networks.keys()]*n_resolutions, 
                            "cluster_resolution":clusterResolutions*n_networks,
                            "num_clusters":num_clusters*n_networks,
                            "p":[1]*n_networks*n_resolutions,
                            "threshold_number":[int(network_sizes['numEdges'].max())]*n_networks*n_resolutions,
                            "pruning":["none"]*n_networks*n_resolutions})
  predictions = {
    i: evaluator.trainCausalModelAndPredict(expression=train_data,
                                  baseNetwork=networks[experiments.loc[i,'network']],
                                  memoizationName=outputs + "/" + str(i) + ".celloracle.oracle", 
                                  perturbations=perturbationsToPredict,
                                  clusterColumnName = experiments.loc[i, "cluster_resolution"],
                                  pruningParameters = {"p":experiments.loc[i,'p'], 
                                                        "threshold_number":experiments.loc[i,'threshold_number']}) 
    for i in experiments.index
  }
  other = None
  return experiments, predictions, other


def plot(evaluationResults, output):
  """Plots specific to this experiment.

  Args:
      evaluationResults (_type_): dataframe with evaluation results, often one row per combo of method and gene
      output (_type_): where to save output
  """