PROJECT_PATH = '/home/gary/cahan_rotation/'
import os
import shutil
import unittest
import pandas as pd
import numpy as np
import scanpy as sc 
import torch
import torch.nn as nn
import anndata
import networkx as nx
import matplotlib.pyplot as plt

import sys
import importlib

os.chdir(PROJECT_PATH + "perturbation_data")
sys.path.append("setup")
import ingestion
importlib.reload(ingestion)

os.chdir(PROJECT_PATH + "perturbation_benchmarking")
sys.path.append("src")
import ggrn
importlib.reload(ggrn)
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import load_networks
importlib.reload(load_networks)

import warnings
warnings.filterwarnings("ignore")

np.set_printoptions(edgeitems=30, linewidth=100000, 
    formatter=dict(float=lambda x: "%6.4f" % x))

def forward(x, n):
    global adjMat, weight
    x = np.matmul(x, adjMat * weight) + n
    return np.squeeze(np.asarray(x))

def simulateKO(control_expression: np.ndarray, 
               noise: np.ndarray, 
               KO_gene_idx: int=0, 
               KO_gene_value: float = 0, 
               maxiter: int=100):
    x = control_expression
    for i in range(100):
        xold = x.copy()
        x[KO_gene_idx] = KO_gene_value
        x = forward(x, noise) 
        if np.linalg.norm(xold - x) < 1e-12:
            break
    x[KO_gene_idx] = KO_gene_value
    return x

G = nx.DiGraph()
G.add_edges_from([(0,1), (1,2), 
                  (2,3), (3,4)])
print(f"G is a DAG? {nx.is_directed_acyclic_graph(G)}")

adjMat = nx.adjacency_matrix(G, range(len(G.nodes))).todense()
adjMat = torch.from_numpy(adjMat)
print("Adjacency Matrix")
print(adjMat.numpy())

weight = adjMat
print("Weight")
print(weight.numpy())

numInst    = 3000
trainShape = (numInst, adjMat.shape[0])
cntrl      = np.zeros(trainShape)

# Create a baseline control expression
source     = np.array([5, 0, 0, 0, 0])
print("\nBias expression level at SS (run simulations to completion 3 times):")
bias = np.ones(5) * 5
bias = simulateKO(bias, source, KO_gene_idx=0, KO_gene_value=source[0])
bias = simulateKO(bias, source, KO_gene_idx=0, KO_gene_value=source[0])
bias = simulateKO(bias, source, KO_gene_idx=0, KO_gene_value=source[0])

# Perturb each gene separately
cntrl += bias
fold=5
for f in range(fold):
    for i in range(numInst//fold*f, numInst//fold*(f+1)):
        cntrl[i] = simulateKO(cntrl[i], source, KO_gene_idx=f, KO_gene_value=np.linspace(1e-2, bias[f], numInst//fold+1)[i-numInst//fold*f])
        
# Add a control expression
cntrl = np.vstack([cntrl, bias])

# Create an anndata and add all the required features
adata = anndata.AnnData(X=cntrl, 
                        dtype=np.float64,
                        obs=pd.DataFrame([f"Sample{i}" for i in range(cntrl.shape[0])], columns=["sample"]),
                        var=pd.DataFrame([f"Gene{i}"   for i in range(cntrl.shape[1])], columns=["gene"]))
adata.obs["perturbation"] = (["Gene0"] * (numInst//fold) + 
                             ["Gene1"] * (numInst//fold) + 
                             ["Gene2"] * (numInst//fold) + 
                             ["Gene3"] * (numInst//fold) +
                             ["Gene4"] * (numInst//fold) +
                             ["Control"])
adata.obs["is_control"]   = adata.obs["perturbation"] == "Control"
adata.obs["is_control_int"]   = [1 if i else 0 for i in adata.obs["is_control"]]
adata.var.index = adata.var.gene
adata.var.index.name = "geneName"
adata.obs.index = [str(i) for i in range(numInst+1)]
adata.obs["consistentW/Perturbation"] = True
adata.obs["logFC"] = -999
adata.obs["spearmanCorr"] = -999

adata.raw = adata

adata.obs["fake_cluster"] = "all_one_cluster"
adata.obs.fake_cluster = adata.obs.fake_cluster.astype("category")
adata.uns["fake_cluster_colors"] = ['#1f77b4']

adata.X = adata.raw.X.copy()

perturbed_and_measured_genes = adata.var.index
perturbed_but_not_measured_genes = list()
print("These genes were perturbed but not measured:")
print(perturbed_but_not_measured_genes)

adata.uns["perturbed_and_measured_genes"] = list(set(adata[~adata.obs.is_control].obs.perturbation))
adata.uns["perturbed_but_not_measured_genes"] = list(perturbed_but_not_measured_genes)

adata = ingestion.describe_perturbation_effect(adata, "knockdown")

adata.uns["weight"] = weight.numpy()
adata.uns["weight_mask"] = adjMat.numpy()
adata.uns["bias"]   = bias
print(adata)

test_data = adata

class TestDCDFG(unittest.TestCase):
    
    os.environ['WANDB_MODE'] = 'offline'
    os.environ["WANDB_SILENT"] = "true"

    def test_NOTEARS_run(self):
        grn = ggrn.GRN(train=test_data, tf_list=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-exp-linear-False",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 20, 
                "num_fine_epochs": 10,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64
                }
        )
        
    def test_NOTEARSLR_run(self):
        grn = ggrn.GRN(train=test_data, tf_list=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-spectral_radius-linearlr-False",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 20, 
                "num_fine_epochs": 10,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64
                }
        )
        
    def test_NOBEARS_run(self):
        grn = ggrn.GRN(train=test_data, tf_list=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-exp-linear-True",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 20, 
                "num_fine_epochs": 10,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64
                }
        )  
        
    def test_DCDFG_run(self):
        grn = ggrn.GRN(train=test_data, tf_list=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-spectral_radius-mlplr-False",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 20, 
                "num_fine_epochs": 10,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64,
                "num_modules": 4,
                }
        )    
    
    def test_NOTEARS_run_and_predict(self):
        grn = ggrn.GRN(train=test_data, tf_list=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-exp-linear-False",
            # method = "DCDFG-spectral_radius-mlplr-False",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 600, 
                "num_fine_epochs": 100,
                "num_gpus": [1] if torch.cuda.is_available() else 0,
                "train_batch_size": 64,
                "regularization_parameter": 0.1,
                }
        )

        posWeight = (grn.models.model.module.weights * grn.models.model.module.weight_mask).detach().cpu().numpy()
        posWeightAnswer = np.array(
            [[0, 1, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [0, 0, 0, 1, 0],
             [0, 0, 0, 0, 1],
             [0, 0, 0, 0, 0]])
        
        bias = grn.models.model.module.biases.detach().cpu().numpy()
        biasAnswer = np.array([5.0, 0.0, 0.0, 0.0, 0.0])

        np.testing.assert_almost_equal(posWeight, posWeightAnswer, decimal=2)
        np.testing.assert_almost_equal(bias, biasAnswer, decimal=2)
    
        control = test_data.X[-1,:]
        koAnswer2 = grn.models.model.simulateKO(control.copy(), 
                                                [0, 1, 2, 3, 4], 
                                                [0.0, 0.0, 0.0, 0.0, 0.0])
        print(koAnswer2)
        
if __name__ == '__main__':
    unittest.main()
