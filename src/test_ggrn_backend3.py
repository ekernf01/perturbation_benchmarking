


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

example_data = sc.read_h5ad("../accessory_data/nakatake.h5ad")
example_network = load_networks.LightNetwork(files=["../accessory_data/human_promoters.parquet"])
example_perturbations = (("KLF8", 0), ("GATA1", 0), ("empty_vector", np.nan))

def simulate_autoregressive(num_controls = 2, num_features = 5, seed = 0, num_steps = 1):
    # This is a pretty normal simulation combining controls and perturbed samples, with linear propagation of perturbation effects.
    # The weirdest thing here is that we ensure the very first control
    # is at steady state, and we tell the algorith so.
    # This helps when S=2 because otherwise the loss has only terms like
    # y = GGx, and GGx = HHx for weird choices of H. Like H=-G or H=P^TGP for a permutation matrix P.
    metadata = pd.DataFrame({
        "matched_control":                      np.tile([i for i in range(num_controls)], num_features + 1),
        "is_control":                           np.repeat([True,  False], [num_controls, num_features*num_controls]),
        "is_treatment":                         np.repeat([False, True ], [num_controls, num_features*num_controls]),
        "is_steady_state":                      [True] + [False]*((num_features+1)*num_controls-1),
        "perturbation":                         np.repeat(["-999"] + [str(i) for i in range(num_features) ], num_controls),
        "expression_level_after_perturbation":  np.repeat(10, (num_features + 1)*num_controls),
    })
    metadata.index = [str(i) for i in metadata.index]
    np.random.seed(seed)
    all_controls = np.random.random((num_controls, num_features))
    G = np.random.random((num_features,num_features))
    # Make sure a steady state exists, then move a control to it via power method
    G = G / np.max(np.real(np.linalg.eigvals(G))) 
    for i in range(100):
        all_controls[0,:] = G.dot(all_controls[0,:])

    def perturb(one_control, perturbation, expression_level_after_perturbation):
        x = one_control.copy()
        if int(perturbation) >= 0:
            x[int(perturbation)] = float(expression_level_after_perturbation)
        return x
    def advance(one_control, perturbation, expression_level_after_perturbation, num_steps):
        one_control = perturb(one_control, perturbation, expression_level_after_perturbation)
        for _ in range(num_steps):
            one_control = G.dot(one_control)
            one_control = perturb(one_control, perturbation, expression_level_after_perturbation)
        return one_control

    linear_autoregressive = anndata.AnnData(
        dtype = np.float32,
        X = np.column_stack(
            [all_controls[i, :] for i in range(num_controls)] + 
            [
                advance(
                    all_controls[metadata.loc[i, "matched_control"], :],
                    metadata.loc[i, "perturbation"], 
                    metadata.loc[i, "expression_level_after_perturbation"], 
                    num_steps = num_steps
                )
                for i in metadata.index if metadata.loc[i, "is_treatment"]
            ]
        ).T,
        obs = metadata,
    )
    return linear_autoregressive, G

class TestBackend3(unittest.TestCase):
    def test_matching_and_loading(self):
        td = ggrn_backend3.ggrn_backend3.MatchControls(example_data, "random")
        # td = ggrn_backend3.ggrn_backend3.MatchControls(train_data, "closest")
        torchy_dataset = ggrn_backend3.ggrn_backend3.AnnDataMatchedControlsDataSet(td, "user")
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

    def test_correctness(self):
        for S in [1,2]:
            linear_autoregressive, G = simulate_autoregressive(num_controls=10, num_features = 3, num_steps=S)
            model = ggrn_backend3.ggrn_backend3.GGRNAutoregressiveModel(linear_autoregressive, matching_method = "user")
            model.train(
                S = S,
                regression_method = "linear",
                low_dimensional_structure = "none",
                low_dimensional_training = None,
                network = None,     
                max_epochs=10000,
                learning_rate=0.001,
                regularization_parameter = 0,
                optimizer = "L-BFGS",
            )  
            [
                np.testing.assert_almost_equal(p.detach().numpy(), G, decimal=1)
                for n,p in model.model.named_parameters() if n=="G.weight"
            ]



    def test_user_interface(self):
        linear_autoregressive, G = simulate_autoregressive(num_controls=10, num_features = 3)
        for matching_method in [
            "user",
            # "random", 
            # "closest",
        ]:
            for regression_method in ["linear"]: #, "multilayer_perceptron"]:
                for low_dimensional_structure in ["none"]: #, "RGQ"]:
                    if low_dimensional_structure == "RGQ":
                        low_dimensional_trainings = ["supervised", "PCA", "fixed"]
                    else:
                        low_dimensional_trainings = [None]
                    for low_dimensional_training in low_dimensional_trainings:
                        model = ggrn_backend3.ggrn_backend3.GGRNAutoregressiveModel(linear_autoregressive, matching_method = matching_method)
                        model.train(
                            S = 1,
                            regression_method = regression_method,
                            low_dimensional_structure = low_dimensional_structure,
                            low_dimensional_training = low_dimensional_training,
                            network = example_network if low_dimensional_training=="fixed" else None,
                            optimizer="L-BFGS",
                        )
                        y = model.predict(example_perturbations)

if __name__ == '__main__':
    unittest.main()
