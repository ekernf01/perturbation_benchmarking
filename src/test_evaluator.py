import unittest
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import shutil

# Deal with various modules specific to this project
import importlib
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking/")
sys.path.append(os.path.expanduser(PROJECT_PATH + 'network_collection/load_networks'))
sys.path.append(os.path.expanduser(PROJECT_PATH + 'perturbation_data/load_perturbations')) 
sys.path.append(os.path.expanduser(PROJECT_PATH + 'perturbation_benchmarking/src')) 
import evaluator
import load_networks
import load_perturbations
importlib.reload(evaluator) 
importlib.reload(load_networks) 
importlib.reload(load_perturbations)
os.environ["GRN_PATH"]           = PROJECT_PATH + "network_collection/networks"
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbation_data/perturbations"
test_expression = sc.read_h5ad(os.environ["PERTURBATION_PATH"] + "/software_test/" + "test.h5ad")
test_network = pd.DataFrame({"peak": np.nan, "gene": ["AATF", "ALX3", "MYOD1"], "AATF": [1.0,0,1], "ALX3": [1.0,1, 0], "MYOD1": [0.0,0,1]})
test_perturbations = pd.DataFrame(
                    [
                        ("AATF", 0.0), # no connections
                        ("ALX3",0.0),  # not enough connections
                        ("MYOD1",0.0)  # enough connections
                ], 
                columns = ["perturbation", "expression_level_after_perturbation"]
            )
test_expression.obs["expression_level_after_perturbation"] = 0

class TestNetworkHandling(unittest.TestCase):

    def test_makeNetworkDense(self):
        pd.testing.assert_frame_equal(
            test_network.copy(), 
            evaluator.makeNetworkDense(
                evaluator.makeNetworkSparse(test_network, 0.0)
            )
        )
        X = evaluator.makeNetworkDense(
            evaluator.makeRandomNetwork(density=1,  TFs=["AATF", "ALX3", "MYOD1"], targetGenes = ["AATF", "ALX3", "MYOD1"])
        )
    
    def test_pivotNetworkWideToLong(self):
        tallified = evaluator.pivotNetworkWideToLong(
                pd.DataFrame({
                    "gene_short_name": ["a", "b", "c"],
                    "peak": "blerghhh",
                    "A": [1, 0, 0],
                    "B": [1, 1, 0],
                    "C": [1, 0, 1],
                })
            )
        tallified.index = range(5)
        pd.testing.assert_frame_equal(
            tallified,
            pd.DataFrame({
                "regulator": ['A', 'B', 'B', 'C', 'C'],
                "target":    ['a', 'a', 'b', 'a', 'c'],
                "weight": 1,
            })
        )

class TestEvaluation(unittest.TestCase):

    def test_evaluateOnePrediction(self):
        self.assertIsNotNone(
            evaluator.evaluateOnePrediction(
                expression =  test_expression,
                predictedExpression = test_expression, 
                baseline = test_expression["Control",:], 
                doPlots=False, 
                outputs="demo",
                experiment_name="test",    
            )
        )
    
    def test_evaluateCausalModel(self):
        os.makedirs("temp", exist_ok=True)
        self.assertIsNotNone(
            evaluator.evaluateCausalModel(
                heldout     = {"demo": test_expression}, 
                predictions = {"demo": test_expression}, 
                baseline    = {"demo": test_expression["Control",:]}, 
                experiments = pd.DataFrame(["demo"], ["demo"], ["demo"]),
                outputs = "temp", 
                factor_varied = "demo",
                default_level = "demo", 
                classifier = None
            )
        )
        shutil.rmtree("temp")
    
if __name__ == '__main__':
    unittest.main()
