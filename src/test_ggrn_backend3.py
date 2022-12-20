PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import sys 
import os
import ggrn_backend3.ggrn_backend3
import unittest
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import load_networks
import scanpy as sc

class TestModelRuns(unittest.TestCase):
    def test_user_interface(self):
        train_data = sc.read_h5ad("../accessory_data/nakatake.h5ad")
        network = load_networks.LightNetwork(files=["../accessory_data/human_promoters.parquet"])
        example_perturbations = (("KLF8", 0), ("GATA1", 0), ("empty_vector", np.nan))
        for matching in "random", "closest":
            for regression_method in "linear", "multilayer_perceptron":
                for low_dimensional_structure in "none", "QiGQ":
                    if low_dimensional_structure == "QiGQ":
                        low_dimensional_trainings = ["supervised", "PCA", "fixed"]
                    else:
                        low_dimensional_trainings = ["this_is_ignored"]
                    for low_dimensional_training in low_dimensional_trainings:
                        model = ggrn_backend3.AutoregressiveModel()
                        model.train(
                            train_data,
                            is_treatment = None, # defaults to ~train_data.obs["is_control"]
                            is_control   = None, # defaults to train_data.obs["is_control"]
                            is_steady_state = None, # defaults to is_control
                            matching = matching,
                            S = 1,
                            regression_method = regression_method,
                            low_dimensional_structure = low_dimensional_structure,
                            low_dimensional_training = low_dimensional_training,
                            network = network,
                        )
                        y = model.predict(example_perturbations)

if __name__ == '__main__':
    unittest.main()
