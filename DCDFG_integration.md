### Background

Say a method "works" if the predictions are better than the mean or median of the training data. Right now, on Nakatake, none of our stuff works. I suspect, but cannot yet formally show, that DCDFG is one of the only causal network inference methods that will work on the type of data we are dealing with. Our current goal is to figure out where things stop working on the road from "DCDFG on Frangieh data" to "Our stuff on Nakatake data". I suspect that any deviation from DCDFG's data and methods will stop working, so I want to start exactly at, or very close, to their figure 5. 

The [DCDFG repo](https://github.com/Genentech/dcdfg) is here. It implements not just DCD-FG, but also a slate of related methods with names like NOTEARS-LR and NOBEARS. Ask Eric for a rundown of the similarities and differences some time. 

### Integrating their code 

- Download the [DCDFG repo](https://github.com/Genentech/dcdfg). Re-run their benchmarks starting from `run_perturbseq_linear.py`. One run of that script produces one dot from one panel of figure 5. If computation time becomes an obstacle: fewer replicates are fine, and the most important task is to show that DCDFG beats NOTEARS-LR.
- The feature to simulate KO's did not exist originally. In `linear_baseline/model.py` and `lowrank_linear_baseline/model.py` and `lowrank_mlp/model.py`, I added the following function.  

        def simulateKO(self, control_expression: np.ndarray, KO_gene_idx: int, KO_gene_value: float = 0, maxiter=100, maxiter_cyclic = 1):
            """Simulate a knockout experiment outcome given control expression and given which gene is KO'd."""
            if not self.module.check_acyclicity():
                print(f"Warning: graph is not acyclic. Predictions may diverge (give NaN's). Setting maxiter to {maxiter_cyclic}.")
                maxiter = maxiter_cyclic
            x = torch.from_numpy(control_expression)
            x = x.float()
            for i in range(maxiter):
                xold = x
                x[:, KO_gene_idx] = KO_gene_value
                x = self.module.forward(x)
                x[:, KO_gene_idx] = KO_gene_value
                if torch.linalg.vector_norm(xold - x) < 1e-12:
                    break
            return x.detach().numpy()

- The above simulator is not well tested. We gotta fix it. It sometimes works, especially for `model="linear"` and `"linearlr"`, but the quality of the results is unclear. For `model="mlplr"`, it returns almost all nans. I believe the nans are due to numerical overflow driven by cycles in the causal graph. Some possible solutions:
    - Make sure the causal graph is acyclic once training is done.
    - Set `maxiter` lower so that even with cycles it won't blow up. (CellOracle has no method to prevent cycles and would set `maxiter` to 3.)
- Pull the latest commits from the [benchmarking repo](https://github.com/ekernf01/perturbation_benchmarking)Run a basic experiment (suggest: "test" and "1.0_1") and make sure it works. Inspect the basic structure of the code and output. Ask Eric questions. 
- There is a copy of DCDFG already mostly integrated with ggrn. What remains to be done:
    - Add unit tests to `test_ggrn.py`. We should fit DCD-FG to an extremely small dataset, make predictions, and ensure the predictions are not NaN. Ideally, they would also test correctness using a small linear + Gaussian simulation. See existing examples in that file.
    - Modify the environment to allow dcdfg to run. Backwards compatibility is not a high priority yet because all our experiments so far 1) have relatively simple reqs and 2) have sucked. We can do whatever works here.
    - Get the unit tests to pass. The DCD-FG integration focuses on these files, so expect errors in:
        - `src/dcdfg_wrapper/dcdfg_wrapper.py`
        - `src/ggrn.py`

