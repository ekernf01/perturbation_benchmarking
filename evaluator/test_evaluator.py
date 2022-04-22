import unittest
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co

# Deal with various modules specific to this project
import importlib
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "benchmarking/")
sys.path.append(os.path.expanduser(PROJECT_PATH + 'networks/load_networks'))
sys.path.append(os.path.expanduser(PROJECT_PATH + 'perturbations/load_perturbations')) 
sys.path.append(os.path.expanduser(PROJECT_PATH + 'benchmarking/evaluator')) 
import evaluator
import load_networks
import load_perturbations
importlib.reload(evaluator) 
importlib.reload(load_networks) 
importlib.reload(load_perturbations)
os.environ["GRN_PATH"]           = PROJECT_PATH + "networks/networks"
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbations/perturbations"
test_expression = sc.read_h5ad(os.environ["PERTURBATION_PATH"] + "/software_test/" + "test.h5ad")
test_network = co.data.load_human_promoter_base_GRN()

class TestEvaluation(unittest.TestCase):

    def test_evaluateCausalModel(self):
        self.assertIsNotNone(
            evaluator.evaluateCausalModel(
                test_expression, 
                {'AATF':np.array(test_expression["AATF",:].X)}, #, "ALX3":np.nan
                baseline = test_expression["Control",:].X, 
                doPlots=False
            )
        )

class TestNetworkHandling(unittest.TestCase):

    def test_makeNetworkDense(self):
        pd.testing.assert_frame_equal(
            test_network.copy(), 
            evaluator.makeNetworkDense(
                evaluator.makeNetworkSparse(test_network, 0.0)
            )
        )
        
        X = evaluator.makeNetworkDense(
            evaluator.makeRandomNetwork(density=1,  TFs=["AATF", "ALX3", "MYOD1"])
        )
        self.assertIsNotNone(
            dict(X.apply(lambda x: x[x>0].index.values, axis=1))
        )

class TestModelTraining(unittest.TestCase):

    def test_trainCausalModelAndPredict(self):
        self.assertIsNotNone(
            evaluator.trainCausalModelAndPredict(
                test_expression, 
                test_network, 
                memoizationName = None,
                perturbations = [("AATF", 0.0), # no connections
                                  ("ALX3",0.0), # not enough connections
                                  ("MYOD1",0.0) # enough connections
                                  ],
                clusterColumnName = "fake_cluster"
            )
        )
    def test_sparseBaseNetwork(self):
        self.assertIsNotNone(
            evaluator.trainCausalModelAndPredict(
                test_expression, 
                evaluator.makeRandomNetwork(density=1, TFs=["AATF", "ALX3", "MYOD1"]), 
                memoizationName = None,
                perturbations = [("AATF", 0.0), 
                                  ("ALX3",0.0), 
                                  ("MYOD1",0.0) 
                                  ],
                clusterColumnName = "fake_cluster"
            )
        )            
            
if __name__ == '__main__':
    unittest.main()


  
    
    
    
    
# def networkEdgesToMatrix(networkEdges, regulatorColumn=0, targetColumn=1):
#     """Reformat a network from a two-column dataframe to the way that celloracle needs its input."""
#     X = pd.crosstab(networkEdges.iloc[:,targetColumn], networkEdges.iloc[:,regulatorColumn])
#     del networkEdges
#     gc.collect()
#     X = 1.0*(X > 0)
#     X = X.rename_axis('gene_short_name').reset_index()
#     X = X.rename_axis('peak_id').reset_index()
#     # CellOracle's preferred format wastes gobs of memory unless you sparsify.
#     X.iloc[:,2:] = X.iloc[:,2:].astype(pd.SparseDtype("float", 0))
#     gc.collect()
#     return X

# def countMatrixEdges(network):
#     """Count the number of connections in a network that is formatted how CO expects it to be formatted."""
#     return 1.0*network.iloc[:,2:].sum().sum()

# humanTFs = pd.read_csv("../accessory_data/humanTFs.csv")
# targetGenes = co.data.load_human_promoter_base_GRN()["gene_short_name"]

# def makeRandomNetwork(density = 0, seed = 0, TFs = humanTFs, targetGenes = targetGenes ):
#     """Generate a random network formatted the way that celloracle needs its input."""
#     np.random.seed(seed)
#     X = pd.DataFrame(
#             np.random.binomial(
#                 n = 1, 
#                 p=density,
#                 size=(
#                     len(targetGenes), 
#                     len(TFs['HGNC symbol'])
#                 )
#             ),
#             columns = TFs['HGNC symbol'], 
#             index = targetGenes
#         )
#     X.rename_axis('gene_short_name', inplace=True)
#     X.reset_index(inplace=True)
#     X.rename_axis('peak_id', inplace=True)
#     X.reset_index(inplace=True)
#     # CellOracle's preferred format wastes gobs of memory unless you sparsify.
#     X.iloc[:,2:] = X.iloc[:,2:].astype(pd.SparseDtype("float", round(density)))
#     gc.collect()
#     return X

# def splitData(adata, allowedRegulators, minTestSetSize = 250, perturbationColName = "perturbation"):
#     """Determine a train-test split satisfying constraints imposed by base networks and available data.
    
# A few factors complicate the training-test split. 

# - Perturbed genes may be absent from most base GRN's due to lack of motif information or ChIP data. These are excluded from the test data to avoid obvious failure cases.
# - Perturbed genes may not be measured. These are excluded from test data because we don't know to what extent they were overexpressed.

# In both cases, we still use those perturbed profiles as training data, hoping they will provide useful info about attainable cell states and downstream causal effects. 

# For some collections of base networks, there are many factors ineligible for use as test data -- so many that we use all the eligible ones for test and the only ineligible ones for training. For other cases, such as dense base networks, we have more flexibility, so we send some perturbations to the training set at random even if we would be able to use them in the test set.

# parameters:
#     - adata: AnnData object satisfying the expectations outlined in the accompanying collection of perturbation data.
#     - allowedRegulators: list or set of features allowed to be in the test set. In CellOracle, this is usually constrained by motif/chip availability. 

#     """
#     allowedRegulators = set(allowedRegulators)
#     allowedRegulators = allowedRegulators.difference(adata.uns["perturbed_but_not_measured_genes"])
#     if len(allowedRegulators) < minTestSetSize:
#         raise ValueError(f"minTestSetSize too big; set to {minTestSetSize} but {len(allowedRegulators)} available.")
#     testSetPerturbations     = set(adata.obs[perturbationColName]).intersection(allowedRegulators)
#     trainingSetPerturbations = set(adata.obs[perturbationColName]).difference(allowedRegulators)
#     if len(trainingSetPerturbations) < minTestSetSize:
#         swap = np.random.default_rng(seed=0).choice(list(testSetPerturbations), 
#                                                minTestSetSize - len(trainingSetPerturbations), 
#                                                replace = False)
#         testSetPerturbations = testSetPerturbations.difference(swap)
#         trainingSetPerturbations = trainingSetPerturbations.union(swap)
#     adata_heldout  = adata[adata.obs[perturbationColName].isin(testSetPerturbations),    :]
#     adata_train    = adata[adata.obs[perturbationColName].isin(trainingSetPerturbations),:]
#     adata_train.obs['perturbation'].unique()
#     perturbationsToPredict = [(gene, adata_heldout[sample, gene].X[0,0]) for sample,gene in 
#                               enumerate(adata_heldout.obs[perturbationColName])] 
#     print("Example perturbations formatted as \n (gene, expression after perturbation)")
#     print(perturbationsToPredict[0:5])
#     print("Test set size:")
#     print(len(testSetPerturbations))
#     print("Training set size:")
#     print(len(trainingSetPerturbations))
#     return adata_train, adata_heldout, perturbationsToPredict
