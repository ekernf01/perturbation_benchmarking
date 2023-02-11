PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import os
import unittest
import pandas as pd
import numpy as np
import scipy
import torch
import anndata
import networkx as nx
import sys
import importlib
import warnings
# Our modules
sys.path.append(os.path.join(PROJECT_PATH, "perturbation_data", "setup"))
import ingestion
importlib.reload(ingestion)
sys.path.append(os.path.join(PROJECT_PATH, "perturbation_benchmarking", "src"))
import ggrn
importlib.reload(ggrn)
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import load_networks
importlib.reload(load_networks)

# Global settings for this test script
warnings.filterwarnings("ignore")
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
np.set_printoptions(edgeitems=30, linewidth=100000, 
    formatter=dict(float=lambda x: "%6.4f" % x))

# Simple simulation code
def forward(x, n):
    global adjMat, weight
    x = np.matmul(x, adjMat * weight) + n
    return np.squeeze(np.asarray(x))

# Unlike the prediction methods we added to DCD-FG, this one iterates all the way to steady state. 
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

# Create a DAG
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
                             ["Control,Gene3"] * (numInst//fold) +
                             ["Gene4,Control"] * (numInst//fold) +
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
adata.X = adata.raw.X.copy()
adata.X = scipy.sparse.csr_matrix(adata.X)
perturbed_and_measured_genes = adata.var.index
perturbed_but_not_measured_genes = list()
print("These genes were perturbed but not measured:")
print(perturbed_but_not_measured_genes)
adata.uns["perturbed_and_measured_genes"] = list(set(adata[~adata.obs.is_control].obs.perturbation))
adata.uns["perturbed_but_not_measured_genes"] = list(perturbed_but_not_measured_genes)
adata = ingestion.describe_perturbation_effect(adata, "knockdown", multiple_genes_hit = True)
adata.uns["weight"] = weight.numpy()
adata.uns["weight_mask"] = adjMat.numpy()
adata.uns["bias"]   = bias
print(adata)
test_data = adata

class TestDCDFG(unittest.TestCase):
    
    os.environ['WANDB_MODE'] = 'offline'
    os.environ["WANDB_SILENT"] = "true"

    def test_NOTEARS_run(self):
        grn = ggrn.GRN(train=test_data, eligible_regulators=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-exp-linear-False",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 2, 
                "num_fine_epochs": 1,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64,
                "verbose": False,
                }
        )
        
    def test_NOTEARSLR_run(self):
        grn = ggrn.GRN(train=test_data, eligible_regulators=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-spectral_radius-linearlr-False",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 20, # This fails with 2 epochs 
                "num_fine_epochs": 1,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64,
                "verbose": False,
                }
        )
        
    def test_NOBEARS_run(self):
        grn = ggrn.GRN(train=test_data, eligible_regulators=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-exp-linear-True",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 2, 
                "num_fine_epochs": 1,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64,
                "verbose": False,
                }
        )  
        
    def test_DCDFG_run(self):
        grn = ggrn.GRN(train=test_data, eligible_regulators=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-spectral_radius-mlplr-False",
            cell_type_sharing_strategy = "identical",
            network_prior = "ignore",
            kwargs = { 
                "num_train_epochs": 2, 
                "num_fine_epochs": 1,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64,
                "verbose": False,
                "num_modules": 4,
            }
        )
    

    def test_NOTEARS_prediction(self):
        grn = ggrn.GRN(train=test_data, eligible_regulators=test_data.var_names, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "DCDFG-spectral_radius-linear-False", 
            kwargs = { 
                "num_train_epochs": 2, 
                "num_fine_epochs": 1,
                "num_gpus": 1 if torch.cuda.is_available() else 0,
                "train_batch_size": 64,
                "verbose": False,
                }
        )

        _ = grn.models.model.module.load_state_dict(
            {
                "weights": torch.nn.Parameter(torch.tensor(np.array(
                [[0.0,   1, 0, 0, 0],
                [0,     0, 1, 0, 0],
                [0,     0, 0, 1, 0],
                [0,     0, 0, 0, 1],
                [0,     0, 0, 0, 0]])), requires_grad=False)
            },  
            strict=False)

        grn.models.model.module.biases = torch.nn.Parameter(torch.tensor(np.array(
            [5.0,    0, 0, 0, 0])), requires_grad=False)
        grn.models.model.module.weight_mask = torch.nn.Parameter(torch.tensor(np.ones((5,5))))
        grn.models.model.module.weights = torch.nn.Parameter(torch.tensor(np.array(
            [[0.0,   1, 0, 0, 0],
             [0,     0, 1, 0, 0],
             [0,     0, 0, 1, 0],
             [0,     0, 0, 0, 1],
             [0,     0, 0, 0, 0]])), requires_grad=False)
        control = test_data.X[-1,:]
        koAnswer2 = grn.predict(
            starting_expression = anndata.AnnData(
                X=np.tile(control.copy(), (7,1)), 
                var = pd.DataFrame(index=[f'Gene{i}' for i in range(5)]),
            ),
            perturbations = [(f'Gene{i}',1000) for i in range(5)] + [("control", 0)] + [("Gene0", np.nan)]
        )
        np.testing.assert_almost_equal(
            koAnswer2.X / np.array([
                [1000, 1000, 5, 5, 5],
                [5, 1000, 1000, 5, 5],
                [5, 5, 1000, 1000, 5],
                [5, 5, 5, 1000, 1000],
                [5, 5, 5,    5, 1000],
                [5, 5, 5,    5,    5],
                [5, 5, 5,    5,    5],
            ]),
            np.ones((7,5)),
            decimal=2
        )

        koAnswer3 = grn.predict(
            perturbations = [(f'Gene{i}',1000) for i in range(5)] + [("control", 0)] + [("Gene0", np.nan)]
        )
        np.testing.assert_almost_equal(
            koAnswer3.X / np.array([
                [1000, 1000, 5, 5, 5],
                [5, 1000, 1000, 5, 5],
                [5, 5, 1000, 1000, 5],
                [5, 5, 5, 1000, 1000],
                [5, 5, 5,    5, 1000],
                [5, 5, 5,    5,    5],
                [5, 5, 5,    5,    5],
            ]),
            np.ones((7,5)),
            decimal=2
        )


    def test_NOTEARS_inference(self):
        grn = ggrn.GRN(train=test_data, eligible_regulators=test_data.var_names, validate_immediately=False)
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
                "verbose": False,
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


if __name__ == '__main__':
    unittest.main()
