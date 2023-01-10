


import os
try:
    PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
    os.chdir(os.path.join(PROJECT_PATH, "perturbation_benchmarking"))
except FileNotFoundError:
    PROJECT_PATH = 'gary_path'
    os.chdir(os.path.join(PROJECT_PATH, "perturbation_benchmarking"))

import sys 
import unittest
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import ggrn_backend3.ggrn_backend3
import load_networks
import scanpy as sc
import numpy as np
import pandas as pd
import torch
import anndata

train_data = sc.read_h5ad("../accessory_data/nakatake.h5ad")
network = load_networks.LightNetwork(files=["../accessory_data/human_promoters.parquet"])
example_perturbations = (("KLF8", 0), ("GATA1", 0), ("empty_vector", np.nan))

def simple_simulation():
    metadata = pd.DataFrame({
        "is_control":      [True]   + [False  for _ in range(20) ],
        "is_treatment":    [False]  + [True   for _ in range(20) ],
        "is_steady_state": [False]  + [False  for _ in range(20) ],
        "perturbation":    ["-999"] + [str(i) for i in range(10) ] + [str(i) for i in range(10) ],
        "expression_level_after_perturbation":   ["0"] + ["0" for _ in range(10) ] + ["10" for _ in range(10) ],
    })

    control = np.random.random((10,1))
    def perturb(control, perturbation, expression_level_after_perturbation):
        x = control.copy()
        x[int(perturbation)] = float(expression_level_after_perturbation)
        return x

    W = np.random.random((10,10))
    linear_one_step = anndata.AnnData(
        dtype = np.float32,
        X = np.column_stack(
            [control] + [
                perturb(
                    W.dot(
                        perturb(
                            control,
                            metadata.loc[i, "perturbation"], 
                            metadata.loc[i, "expression_level_after_perturbation"], 
                        )
                    ),
                    metadata.loc[i, "perturbation"], 
                    metadata.loc[i, "expression_level_after_perturbation"], 
                )
                for i in metadata.index if i != 0
            ]
        ).T,
        obs = metadata,
    )
    return linear_one_step, W

class TestBackend3(unittest.TestCase):
    def test_matching_and_loading(self):
        td = ggrn_backend3.ggrn_backend3.MatchControls(train_data, "random")
        # td = ggrn_backend3.ggrn_backend3.MatchControls(train_data, "closest")
        torchy_dataset = ggrn_backend3.ggrn_backend3.AnnDataMatchedControlsDataSet(td, "random")
        torchy_dataloader = torch.utils.data.DataLoader(torchy_dataset)
        x = torchy_dataset[0]
        batch = next(iter(torchy_dataloader))
        assert "treatment" in x.keys()
        assert "matched_control" in x.keys()
        assert "expression" in x["treatment"].keys()
        assert "metadata"   in x["treatment"].keys()
        assert "expression" in x["matched_control"].keys()
        assert "metadata"   in x["matched_control"].keys()
        assert all(
            [f in x["treatment"]["metadata"].keys() 
            for f in {"is_control", "is_treatment", "is_steady_state", "perturbation_index", "expression_level_after_perturbation"}]
        )

    def test_user_interface(self):
        for matching in [
            "random", 
            # "closest",
        ]:
            for regression_method in ["linear"]: #, "multilayer_perceptron"]:
                for low_dimensional_structure in ["none"]: #, "RGQ"]:
                    if low_dimensional_structure == "RGQ":
                        low_dimensional_trainings = ["supervised", "PCA", "fixed"]
                    else:
                        low_dimensional_trainings = [None]
                    for low_dimensional_training in low_dimensional_trainings:
                        model = ggrn_backend3.ggrn_backend3.GGRNAutoregressiveModel(train_data, matching_method = "random")
                        model.train(
                            S = 1,
                            regression_method = regression_method,
                            low_dimensional_structure = low_dimensional_structure,
                            low_dimensional_training = low_dimensional_training,
                            network = network,
                        )
                        y = model.predict(example_perturbations)

    def test_correctness(self):
        linear_one_step, W = simple_simulation()
        model = ggrn_backend3.ggrn_backend3.GGRNAutoregressiveModel(linear_one_step, matching_method = "random")
        model.train(
            S = 1,
            regression_method = "linear",
            low_dimensional_structure = "none",
            low_dimensional_training = None,
            network = None,     
            max_epochs=10000,
            regularization_parameter = 0,
        )
        model
        x = linear_one_step.X[0,:].copy()
        x[0] = 10
        np.dot(W, x) - linear_one_step.X[11,:]
        [p.detach().numpy()-W for n,p in model.model.named_parameters() if n=="G.weight"]






if __name__ == '__main__':
    unittest.main()
