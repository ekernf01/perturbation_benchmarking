PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import os
import unittest
import pandas as pd
import scanpy as sc 
import anndata
os.chdir(PROJECT_PATH + "perturbation_benchmarking")

import sys
import importlib
sys.path.append("src")
import predict
importlib.reload(predict)

train   = sc.read_h5ad("../accessory_data/nakatake.h5ad")
network = pd.read_parquet("../accessory_data/human_promoters.parquet")

class TestModelInstantiation(unittest.TestCase):
    def test_make_GRN(self):
        self.assertIsInstance(
            predict.GRN(train, network = network), predict.GRN
        )
    def test_validate_input(self):
        self.assertTrue(
            predict.GRN(train, network = network).check_perturbation_dataset() 
        )    
    def test_extract_features(self):
        grn    = predict.GRN(train, network = network)
        self.assertIsNone(
            grn.extract_features(method = "tf_rna")
        )    
    def test_simple_fit_and_predict(self):
        grn    = predict.GRN(train[0:100, 0:100].copy(), network = network, validate_immediately=False)
        grn.extract_features(method = "tf_rna")
        self.assertIsNone(
            grn.fit(
                method = "linear", 
                confounders = ["donor", "batch"],
                cell_type_sharing_strategy = "identical", 
                network_prior = "ignore",
                pruning_strategy = "none", 
                pruning_parameter = None, 
                projection = "factor_graph", 
            )
        )
        p = grn.predict((("KLF8", 0), ("GATA1", 0)))
        print(p)
        self.assertIsInstance(
            p, anndata.AnnData
        )

    def test_network_fit_and_predict(self):
        grn    = predict.GRN(train[0:100, 0:100].copy(), network = network, validate_immediately=False)
        grn.extract_features(method = "tf_rna")
        self.assertIsNone(
            grn.fit(
                method = "linear", 
                confounders = ["donor", "batch"], 
                cell_type_sharing_strategy = "identical", 
                network_prior = "restrictive", 
                pruning_strategy = "none", 
                pruning_parameter = None, 
                projection = "factor_graph",     
            )
        )
        p = grn.predict((("KLF8", 0), ("GATA1", 0)))
        print(p)
        self.assertIsInstance(
            p, anndata.AnnData
        )

    def test_pruned_fit_and_predict(self):
        grn    = predict.GRN(train[0:100, 0:100].copy(), network=network, validate_immediately=False)
        grn.extract_features(method = "tf_rna")
        grn.fit(
            method = "linear", 
            confounders = ["donor", "batch"], 
            cell_type_sharing_strategy = "identical", 
            network_prior = "ignore", 
            pruning_strategy = "prune_and_refit",
            pruning_parameter = 5000, 
            projection = "factor_graph",   
        )
        p = grn.predict((("KLF8", 0), ("GATA1", 0)))
        print(p)
        self.assertIsInstance(
            p, anndata.AnnData
        )

if __name__ == '__main__':
    unittest.main()
