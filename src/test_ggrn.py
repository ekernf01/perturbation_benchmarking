PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import os
import shutil
import unittest
import pandas as pd
import numpy as np
import scanpy as sc 
import scipy
import anndata
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import sys
import importlib
sys.path.append("src")
import ggrn
importlib.reload(ggrn)
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import load_networks
importlib.reload(load_networks)



train   = sc.read_h5ad("../accessory_data/nakatake.h5ad")
train.X = scipy.sparse.csr_matrix(train.X) 
network = load_networks.LightNetwork(files=["../accessory_data/human_promoters.parquet"])
example_perturbations = [("Control,KLF8", "0,0"), ("GATA1", 0), ("GATA1,GFP", f"0,{str(np.nan)}"), ("empty_vector", np.nan)]

class TestModelRuns(unittest.TestCase):
    def test_make_GRN(self):
        self.assertIsInstance(
            ggrn.GRN(train, network = network), ggrn.GRN
        )
        self.assertIsInstance(
            ggrn.GRN(train), ggrn.GRN
        )
        grn = ggrn.GRN(train, network = network)

    def test_validate_input(self):
        self.assertTrue(
            ggrn.GRN(train, network = network).check_perturbation_dataset() 
        )    
        self.assertTrue(
            ggrn.GRN(train).check_perturbation_dataset() 
        )    
    def test_extract_tf_activity(self):
        grn    = ggrn.GRN(train, network = network)
        self.assertIsNone(
            grn.extract_tf_activity(method = "tf_rna")
        )    
        grn    = ggrn.GRN(train)
        self.assertIsNone(
            grn.extract_tf_activity(method = "tf_rna")
        )    

    def test_fit_various_methods(self):
        grn    = ggrn.GRN(train[0:100, 0:10].copy(), validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        for regression_method in [
                "mean", 
                "median", 
                "GradientBoostingRegressor",
                "ExtraTreesRegressor",
                "KernelRidge",
                "RidgeCV", 
                "RidgeCVExtraPenalty",
                "LassoCV",
                "ElasticNetCV",
                "LarsCV",
                "LassoLarsIC",
                "OrthogonalMatchingPursuitCV",
                "ARDRegression",
                "BayesianRidge",
                # "BernoulliRBM",
                # "MLPRegressor",
            ]:
            grn.fit(
                method = regression_method, 
                cell_type_sharing_strategy = "identical", 
                network_prior = "ignore",
                pruning_strategy = "none", 
                pruning_parameter = None, 
                do_parallel=False, #catch bugs easier
            )
            x1 = grn.predict(perturbations = example_perturbations)  
            x2 = grn.predict(perturbations = example_perturbations, starting_expression=train[0:len(example_perturbations), 0:10].copy())
    
    def test_simple_fit_and_predict(self):
        grn    = ggrn.GRN(train[0:100, 0:100].copy(), validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore",
            pruning_strategy = "none", 
            pruning_parameter = None, 
        )
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore",
            pruning_strategy = "none", 
            pruning_parameter = None, 
            do_parallel=False,
        )
        shutil.rmtree("tempfolder23587623", ignore_errors=True)
        self.assertTrue(grn.save_models("tempfolder23587623") is None)
        shutil.rmtree("tempfolder23587623")
        p = grn.predict(example_perturbations)
        self.assertIsInstance(
            p, anndata.AnnData
        )

    def test_network_fit_and_predict(self):
        grn    = ggrn.GRN(train[0:100, 0:100].copy(), network = network, validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive", 
            pruning_strategy = "none", 
            pruning_parameter = None, 
        )
        p = grn.predict(example_perturbations)
        self.assertIsInstance(
            p, anndata.AnnData
        )

    def test_simulation_works(self):
        grn    = ggrn.GRN(train[0:100, 0:100].copy(), network = network, validate_immediately=False)
        # network-based
        p = grn.simulate_data(example_perturbations, effects = "uniform_on_provided_network", noise_sd=1)
        self.assertTrue(ggrn.GRN(p, validate_immediately=True))
        # trained model-based
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive", 
            pruning_strategy = "none", 
            pruning_parameter = None, 
        )
        p1 = grn.simulate_data(example_perturbations, effects = "fitted_models")    
        p2 = grn.simulate_data(example_perturbations, effects = "fitted_models")     
        p3 = grn.simulate_data(example_perturbations, effects = "fitted_models", noise_sd=1)     
        np.testing.assert_almost_equal(p1.X, p2.X)
        self.assertTrue(ggrn.GRN(p1, validate_immediately=True))
        self.assertTrue(ggrn.GRN(p2, validate_immediately=True))
        self.assertTrue(ggrn.GRN(p3, validate_immediately=True))

    def test_pruned_fit_and_predict(self):
        grn    = ggrn.GRN(train[0:100, 0:100].copy(), validate_immediately=False)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore", 
            pruning_strategy = "prune_and_refit",
            pruning_parameter = 5000, 
        )
        p = grn.predict(example_perturbations, do_parallel=False)
        self.assertIsInstance(
            p, anndata.AnnData
        )



# Make some simulated data that ridge regressions should be able to get exactly right
easy_simulated = anndata.AnnData(
    dtype=np.float32,
    obs = pd.DataFrame(
        {"is_control": True,
        "perturbation": "control",
        "perturbation_type": "overexpression",
        "cell_type": np.resize([1,2], 1000)},
        index = [str(i) for i in range(1000)],
    ),
    var = pd.DataFrame(
        index = ["g" + str(i) for i in range(5)],
    ),
    X = np.random.rand(1000, 5) - 0.5
)
# One gene is a linear function
easy_simulated.X[:,0] = easy_simulated.X[:,1]*2 
# Another is a cell-type-specific linear function
easy_simulated.X[easy_simulated.obs["cell_type"]==1,4] = easy_simulated.X[easy_simulated.obs["cell_type"]==1,2]*2
easy_simulated.X[easy_simulated.obs["cell_type"]==2,4] = easy_simulated.X[easy_simulated.obs["cell_type"]==2,3]*2
# empty network should fit worse than correct network
# correct network should fit roughly same as dense network
empty_network   = load_networks.LightNetwork(df = pd.DataFrame(index=[], columns=["regulator", "target", "weight"]))
correct_network = load_networks.LightNetwork(df = pd.DataFrame(
    {
        "regulator":["g" + str(i) for i in [1, 2, 3]], 
        "target":   ["g" + str(i) for i in [0, 4, 4]], 
        "weight":1,
    }
    ))
dense_network = load_networks.LightNetwork(df = pd.DataFrame({
    "regulator": ["g" + str(i) for _ in range(5) for i in range(5)],
    "target":    ["g" + str(i) for i in range(5) for _ in range(5)],
    "weight": 1,
    }))

class TestModelExactlyRightOnEasySimulation(unittest.TestCase):
    def test_simple_dense(self):
        grn = ggrn.GRN(easy_simulated, validate_immediately=False, eligible_regulators=easy_simulated.var_names)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore", 
            pruning_strategy = "none",
        )
        p = grn.predict([("g1", 0), ("g1", 1), ("g1", 2)], do_parallel=False)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_simple_dense_parallel_predict(self):
        grn = ggrn.GRN(easy_simulated, validate_immediately=False, eligible_regulators=easy_simulated.var_names)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore", 
            pruning_strategy = "none",
        )
        p = grn.predict([("g1", 0), ("g1", 1), ("g1", 2)], do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_network_dense(self):
        grn = ggrn.GRN(easy_simulated, network = dense_network, validate_immediately=False, eligible_regulators=easy_simulated.var_names)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive",
            pruning_strategy = "none"
        )
        p = grn.predict([("g1", 0), ("g1", 1), ("g1", 2)], do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_network_correct(self):
        grn = ggrn.GRN(easy_simulated, network = correct_network, validate_immediately=False, eligible_regulators=easy_simulated.var_names)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive",
            pruning_strategy = "none",
        )
        p = grn.predict([("g1", 0), ("g1", 1), ("g1", 2)], do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_network_empty(self):
        grn = ggrn.GRN(easy_simulated, network = empty_network, validate_immediately=False, eligible_regulators=easy_simulated.var_names)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive",
            pruning_strategy = "none",
        )
        p = grn.predict([("g1", 0), ("g1", 1), ("g1", 2)], do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,0,0], decimal=1)

    def test_ct_specific(self):
        grn = ggrn.GRN(easy_simulated, validate_immediately=False, eligible_regulators=easy_simulated.var_names)
        grn.extract_tf_activity(method = "tf_rna")
        grn.fit(
            method = "RidgeCVExtraPenalty", 
            cell_type_sharing_strategy = "distinct", 
            cell_type_labels='cell_type',
            network_prior = "ignore",
            pruning_strategy = "none",
        )
        
        for dp in [False, True]:
            # Cell type 1: only gene 2 affects gene 4
            dummy_control = easy_simulated[0:3,:].copy()
            dummy_control.X = dummy_control.X*0
            dummy_control.obs["cell_type"] = 1
            p = grn.predict([("g2", 0), ("g2", 1), ("g2", 2)], starting_expression=dummy_control.copy(), do_parallel=dp)
            np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,2,4], decimal=1)
            p = grn.predict([("g3", 0), ("g3", 1), ("g3", 2)], starting_expression=dummy_control.copy(), do_parallel=dp)
            np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,0,0], decimal=1)
            # Cell type 2: only gene 3 affects gene 4
            dummy_control = easy_simulated[0:3,:].copy()
            dummy_control.X = dummy_control.X*0
            dummy_control.obs["cell_type"] = 2
            p = grn.predict([("g2", 0), ("g2", 1), ("g2", 2)], starting_expression=dummy_control.copy(), do_parallel=dp)
            np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,0,0], decimal=1)
            p = grn.predict([("g3", 0), ("g3", 1), ("g3", 2)], starting_expression=dummy_control.copy(), do_parallel=dp)
            np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,2,4], decimal=1)


if __name__ == '__main__':
    unittest.main()
