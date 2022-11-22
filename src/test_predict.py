PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import os
import unittest
import pandas as pd
import numpy as np
import scanpy as sc 
import anndata
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import sys
import importlib
sys.path.append("src")
import predict
importlib.reload(predict)
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import load_networks
importlib.reload(load_networks)

train   = sc.read_h5ad("../accessory_data/nakatake.h5ad")
network = load_networks.LightNetwork(files=["../accessory_data/human_promoters.parquet"])
class TestModelRuns(unittest.TestCase):
    def test_make_GRN(self):
        self.assertIsInstance(
            predict.GRN(train, network = network), predict.GRN
        )
        self.assertIsInstance(
            predict.GRN(train), predict.GRN
        )
        grn = predict.GRN(train, network = network)

    def test_validate_input(self):
        self.assertTrue(
            predict.GRN(train, network = network).check_perturbation_dataset() 
        )    
        self.assertTrue(
            predict.GRN(train).check_perturbation_dataset() 
        )    
    def test_extract_features(self):
        grn    = predict.GRN(train, network = network)
        self.assertIsNone(
            grn.extract_features(method = "tf_rna")
        )    
        grn    = predict.GRN(train)
        self.assertIsNone(
            grn.extract_features(method = "tf_rna")
        )    

    def test_simple_fit_and_predict(self):
        grn    = predict.GRN(train[0:100, 0:100].copy(), validate_immediately=False)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore",
            pruning_strategy = "none", 
            pruning_parameter = None, 
        )
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore",
            pruning_strategy = "none", 
            pruning_parameter = None, 
            do_parallel=False,
        )
        p = grn.predict((("KLF8", 0), ("GATA1", 0)))
        self.assertIsInstance(
            p, anndata.AnnData
        )

    def test_network_fit_and_predict(self):
        grn    = predict.GRN(train[0:100, 0:100].copy(), network = network, validate_immediately=False)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive", 
            pruning_strategy = "none", 
            pruning_parameter = None, 
        )
        p = grn.predict((("KLF8", 0), ("GATA1", 0)))
        self.assertIsInstance(
            p, anndata.AnnData
        )
        p = grn.simulate_data((("KLF8", 0), ("GATA1", 0)), effects = "uniform_on_provided_network", noise_sd=1)
        self.assertIsInstance(
            p, anndata.AnnData
        )
        p = grn.simulate_data((("KLF8", 0), ("GATA1", 0)), effects = "fitted_models", noise_sd=1)
        self.assertIsInstance(
            p, anndata.AnnData
        )

    def test_pruned_fit_and_predict(self):
        grn    = predict.GRN(train[0:100, 0:100].copy(), network=network, validate_immediately=False)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore", 
            pruning_strategy = "prune_and_refit",
            pruning_parameter = 5000, 
        )
        p = grn.predict((("KLF8", 0), ("GATA1", 0)))
        self.assertIsInstance(
            p, anndata.AnnData
        )



# Make some simulated data that ridge regressions should be able to get exactly right
easy_simulated = anndata.AnnData(
    obs = pd.DataFrame(
        {"is_control": True,
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
        grn = predict.GRN(easy_simulated, validate_immediately=False, tf_list=easy_simulated.var_names)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore", 
            pruning_strategy = "none",
        )
        p = grn.predict((("g1", 0), ("g1", 1), ("g1", 2)), do_parallel=False)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_simple_dense_parallel_predict(self):
        grn = predict.GRN(easy_simulated, validate_immediately=False, tf_list=easy_simulated.var_names)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore", 
            pruning_strategy = "none",
        )
        p = grn.predict((("g1", 0), ("g1", 1), ("g1", 2)), do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_network_dense(self):
        grn = predict.GRN(easy_simulated, network = dense_network, validate_immediately=False, tf_list=easy_simulated.var_names)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive",
            pruning_strategy = "none"
        )
        p = grn.predict((("g1", 0), ("g1", 1), ("g1", 2)), do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_network_correct(self):
        grn = predict.GRN(easy_simulated, network = correct_network, validate_immediately=False, tf_list=easy_simulated.var_names)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive",
            pruning_strategy = "none",
        )
        p = grn.predict((("g1", 0), ("g1", 1), ("g1", 2)), do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,2,4], decimal=1)

    def test_network_empty(self):
        grn = predict.GRN(easy_simulated, network = empty_network, validate_immediately=False, tf_list=easy_simulated.var_names)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "identical", 
            network_prior = "restrictive",
            pruning_strategy = "none",
        )
        p = grn.predict((("g1", 0), ("g1", 1), ("g1", 2)), do_parallel=True)
        np.testing.assert_almost_equal(p[:, "g0"].X.squeeze(), [0,0,0], decimal=1)

    def test_ct_specific(self):
        grn = predict.GRN(easy_simulated, validate_immediately=False, tf_list=easy_simulated.var_names)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            cell_type_sharing_strategy = "distinct", 
            cell_type_labels='cell_type',
            network_prior = "ignore",
            pruning_strategy = "none",
        )
        # Cell type 1: only gene 2 affects gene 4
        p = grn.predict((("g2", 0), ("g2", 1), ("g2", 2)), starting_states=easy_simulated.obs["cell_type"]==1, do_parallel=False)
        np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,2,4], decimal=1)
        p = grn.predict((("g3", 0), ("g3", 1), ("g3", 2)), starting_states=easy_simulated.obs["cell_type"]==1)
        np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,0,0], decimal=1)
        # Cell type 2: only gene 3 affects gene 4
        p = grn.predict((("g2", 0), ("g2", 1), ("g2", 2)), starting_states=easy_simulated.obs["cell_type"]==2)
        np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,0,0], decimal=1)
        p = grn.predict((("g3", 0), ("g3", 1), ("g3", 2)), starting_states=easy_simulated.obs["cell_type"]==2)
        np.testing.assert_almost_equal(p[:,"g4"].X.squeeze(), [0,2,4], decimal=1)

if __name__ == '__main__':
    unittest.main()
