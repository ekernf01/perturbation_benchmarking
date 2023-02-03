from multiprocessing.sharedctypes import Value
import anndata
from joblib import Parallel, delayed, cpu_count, dump
import pandas as pd
import os
import sklearn.linear_model
import sklearn.ensemble
import sklearn.neural_network
import sklearn.kernel_ridge
import sklearn.dummy
import numpy as np
import gc 
import psutil
import dcdfg_wrapper.dcdfg_wrapper as dcdfg_wrapper 

# Project-specific paths
import sys 
import importlib
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
DEFAULT_TF_LIST = pd.read_csv("../accessory_data/humanTFs.csv")
DEFAULT_TF_LIST = [g for g in DEFAULT_TF_LIST.loc[DEFAULT_TF_LIST["Is TF?"]=="Yes","HGNC symbol"]]
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
import load_networks
importlib.reload(load_networks)
import load_perturbations
importlib.reload(load_perturbations)


class LinearSimulator:
    """Stand-in for sklearn.linear_model.RidgeCV to allow simulation prior to model fitting."""
    def __init__(self, dimension) -> None:
        self.coef_ = np.ones(dimension) / dimension
    def predict(self, X):
        return X.dot(self.coef_)

class GRN:
    """
    Flexible inference of gene regulatory network models.
    """
    
    def __init__(
        self, train: anndata.AnnData, 
        network: load_networks.LightNetwork = None, 
        tf_list = DEFAULT_TF_LIST, 
        validate_immediately = True
        ):
        """Create a GRN object.

        Args:
            train (anndata.AnnData): Training data. Should conform to the requirements in load_perturbations.check_perturbation_dataset().
            network (pd.DataFrame, optional): LightNetwork object containing prior knowledge about regulators and targets.
            tf_list (_type_, optional): List of gene names that are allowed to be regulators.
        """
        self.train = train 
        assert network is None or type(network)==load_networks.LightNetwork
        self.network = network 
        self.tf_list = [tf for tf in tf_list if tf in train.var_names]
        self.models = [None for _ in self.train.var_names]
        self.training_args = {}
        if validate_immediately:
            assert self.check_perturbation_dataset()

    def check_perturbation_dataset(self):
        return load_perturbations.check_perturbation_dataset(ad=self.train)

    def extract_tf_activity(self, train: anndata.AnnData = None, in_place: bool = True, method = "tf_rna"):
        """Create a feature matrix where each row matches a row in self.train 
        each column represents activity of the corresponding TF in self.tf_list.

        Args:
            train: (anndata.AnnData, optional). Expression data to use in feature extraction. Defaults to self.train.
            in_place (bool, optional): If True (default), put results in self.features. Otherwise, return results as a matrix. 
            method (str, optional): How to extract features. Defaults to "tf_rna", which uses the mRNA level of the TF.
        """
        if train is None:
            train = self.train # shallow copy is best here

        if method == "tf_rna":
            print(self.tf_list)
            features = train[:,self.tf_list].X
        else:
            raise NotImplementedError("Only method='tf_rna' is available so far.")
        
        if in_place:
            self.features = features # shallow copy is best here too
        else:
            return features

    def apply_supervised_ml(
        self, 
        FUN, 
        network_prior: str, 
        pruning_strategy: str, 
        pruning_parameter: int,
        do_parallel: bool,
        verbose:int = 1  ,          
    ):
        """Apply a supervised ML method to predict target expression from TF activity.

        Args:
            FUN (function): should accept inputs X,y and return a fitted model, sklearn-style, with a .predict method.
            verbose (int): passes to joblib.Parallel. https://joblib.readthedocs.io/en/latest/parallel.html
            network_prior (str): See args for 'fit' method
            pruning_strategy (str): See args for 'fit' method
            pruning_parameter (int): See args for 'fit' method
            do_parallel (bool): If True, use joblib parallelization to fit models. 
            verbose (int, optional): Passed to joblib.Parallel. Defaults to 1.

        Returns:
            No return value. Instead, this modifies self.models.
        """
        if self.training_args["time_strategy"] == "two_step":
            # Set intermediate expression equal to:
            # - perturbed (for all TFs)
            # - control (for non-TFs)
            intermediate_expression = self.train.X.copy()
            for i,g in enumerate(self.tf_list):
                baseline_expression = self.train[self.train.obs["is_control"],g].X.mean()
                intermediate_expression[:, i] = baseline_expression
            # Set baseline TF activity equal to:
            # - perturbed sample (for just the perturbed TF)
            # - control sample (otherwise)
            baseline_features  = self.features[self.train.obs["is_control"],:].mean(axis=0), 
            baseline_features  = np.repeat(baseline_features, self.train.X.shape[0], axis=0)
            for i,idx in enumerate(self.train.obs.index):
                perturbed_gene = self.train.obs.loc[idx, "perturbation"] 
                if perturbed_gene in self.train.uns["perturbed_and_measured_genes"] and perturbed_gene in self.tf_list:
                    perturbed_gene_index = [g==perturbed_gene for g in self.tf_list]
                    baseline_features[i, perturbed_gene_index] = self.features[i, perturbed_gene_index]
            # Augment data: 
            # - Baseline TF activity predicts intermediate expression 
            # - Perturbed TF activity predicts final expression
            # - Each sample repped twice so metadata also 2x samples
            features_augmented = np.concatenate([baseline_features, self.features])
            target_augmented   = np.concatenate([intermediate_expression,  self.train.X])            
            obs_augmented      = pd.concat([self.train.obs,self.train.obs])
        elif self.training_args["time_strategy"] == "steady_state":
            features_augmented = self.features
            target_augmented   = self.train.X
            obs_augmented      = self.train.obs
        else:
            raise ValueError("'time_strategy' must be 'steady_state' or 'two_step'.")
        
        if do_parallel:     
            def do_loop(network = self.network):
                    m = Parallel(n_jobs=cpu_count()-1, verbose = verbose, backend="loky")(
                        delayed(apply_supervised_ml_one_gene)(
                            train_obs = obs_augmented,
                            target_expr = target_augmented[:,i],
                            features = features_augmented,
                            network = network,
                            training_args = self.training_args,
                            tf_list = self.tf_list, 
                            FUN=FUN, 
                            network_prior=network_prior,
                            target = self.train.var_names[i],
                        )
                        for i in range(len(self.train.var_names))
                    )
                    return m
        else:
            def do_loop(network = self.network):
                return [
                    apply_supervised_ml_one_gene(
                        train_obs = obs_augmented,
                        target_expr = target_augmented[:,i],
                        features = features_augmented,
                        network = network,
                        training_args = self.training_args,
                        tf_list = self.tf_list, 
                        FUN=FUN, 
                        network_prior=network_prior,
                        target = self.train.var_names[i],
                    )
                    for i in range(len(self.train.var_names))
                ]

        if pruning_strategy == "lasso":
            raise NotImplementedError("lasso pruning not implemented yet.")
        elif pruning_strategy == "prune_and_refit":
            print("Fitting")
            self.models = do_loop()
            print("Pruning")
            chunks = []
            for i in range(len(self.train.var_names)):
                if self.training_args["cell_type_sharing_strategy"] == "distinct":
                    for cell_type in self.train.obs[self.training_args["cell_type_labels"]].unique():
                        chunks.append(
                            pd.DataFrame(
                                    {
                                        "regulator": get_regulators(
                                                tf_list = self.tf_list, 
                                                predict_self = self.training_args["predict_self"], 
                                                network_prior = network_prior,
                                                target = self.train.var_names[i], 
                                                network = self.network,
                                                cell_type = cell_type,
                                                ), 
                                        "target": self.train.var_names[i], 
                                        "weight": self.models[i][cell_type].coef_.squeeze(),
                                        "cell_type": cell_type,
                                    } 
                                )
                            )
                elif self.training_args["cell_type_sharing_strategy"] == "identical":
                    chunks.append(
                        pd.DataFrame(
                                {
                                    "regulator": get_regulators(
                                            tf_list = self.tf_list, 
                                            predict_self = self.training_args["predict_self"], 
                                            network_prior = network_prior,
                                            target = self.train.var_names[i], 
                                            network = self.network,
                                            ), 
                                    "target": self.train.var_names[i], 
                                    "weight": self.models[i].coef_.squeeze(),
                                } 
                            )
                        )
                else:
                    raise NotImplementedError("Invalid value of 'cell_type_sharing_strategy' ")

            pruned_network = pd.concat(chunks)
            pruned_network = pruned_network.query("regulator != 'NO_REGULATORS'")
            pruned_network.loc[:, "abs_weight"] = np.abs(pruned_network["weight"])
            pruned_network = pruned_network.nlargest(pruning_parameter, "abs_weight")
            del pruned_network["abs_weight"]
            self.network = load_networks.LightNetwork(df=pruned_network)
            print("Re-fitting")
            self.models = do_loop()                
        elif pruning_strategy is None or pruning_strategy == "none" or np.isnan(pruning_strategy):
            print("Fitting")
            self.models = do_loop() 
        else:
            raise NotImplementedError(f"pruning_strategy should be one of 'none', 'lasso', 'prune_and_refit'; got {pruning_strategy} ")
    
    def simulate_data(
        self,
        perturbations,
        effects: str = "fitted_models",
        noise_sd = None,
        feature_extraction_method = "tf_rna",
        seed = 0,
    ) -> anndata.AnnData:
        """Generate simulated expression data

        Args:
            perturbations (iterable): See GRN.predict()
            effects (str, optional): Either "fitted_models" (use effect sizes from an already-trained GRN) or
                "uniform_on_provided_network" (use a provided network structure with small positive effects for all regulators).
            noise_sd (float, optional): Standard deviation of noise. Defaults to None: noise sd will be extracted from fitted models.
            seed (int): RNG seed. 

        Returns:
            anndata.AnnData: simulated gene expression values.
        """
        if effects == "fitted_models":
            adata = self.predict(perturbations, add_noise=True, noise_sd=noise_sd)
        elif effects == "uniform_on_provided_network":
            if self.network is None:
                raise ValueError("For network-based simulation, network structure must be provided during GRN initialization.")
            self.training_args["predict_self"] = False
            self.training_args["network_prior"] = "restrictive"
            self.training_args["cell_type_sharing_strategy"] = "identical"
            self.training_args["confounders"] = []
            self.training_args["time_strategy"] = "steady_state"
            self.training_args["method"] = "RidgeCV"
            for i in range(len(self.train.var_names)):
                self.models[i] = LinearSimulator(dimension=len(get_regulators(
                    self.tf_list, 
                    predict_self = False, 
                    target = self.train.var_names[i],
                    network = self.network,
                    network_prior=self.training_args["network_prior"],
                )))
            self.extract_tf_activity(method = feature_extraction_method)
            adata = self.predict(perturbations, add_noise=True, noise_sd=noise_sd, seed = seed)
        else:
            raise ValueError("'effects' must be one of 'fitted_models' or 'uniform_on_provided_network'.")
        # Make it pass the usual validation checks, same as all our perturbation data
        adata.obs["perturbation_type"] = self.train.obs["perturbation_type"][0]
        adata.obs["is_control"] = [np.isnan(p[1]) for p in perturbations]
        adata.obs["spearmanCorr"] = np.nan #could simulate actual replicates later if needed
        adata.uns["perturbed_and_measured_genes"]     = [p[0] for p in perturbations if p[0] in adata.var_names]
        adata.uns["perturbed_but_not_measured_genes"] = [p[0] for p in perturbations if p[0] not in adata.var_names]
        adata.raw = adata.copy()
        return adata

    def save_models(self, folder_name:str):
        """Save all the regression models to a folder via joblib.dump (one file per gene).

        Args:
            folder_name (str): Where to save files.
        """
        if ("model" in self.models.__dict__.keys() and
            self.models.model.__class__.__name__ in ["LinearGaussianModel", 
                                                     "LinearGaussianModel_poly", 
                                                     "LinearModuleGaussianModel", 
                                                     "MLPModuleGaussianModel"]):
            raise NotImplementedError(f"Parameter saving/loading is not supported for DCDFG")
        
        os.makedirs(folder_name, exist_ok=True)
        for i,target in enumerate(self.train.var_names):
            dump(self.models[i], os.path.join(folder_name, f'{target}.joblib'))
        return

    def fit(
        self,
        method: str, 
        confounders: list = [], 
        cell_type_sharing_strategy: str = "distinct",   
        cell_type_labels: str = None,
        network_prior: str = "restrictive",    
        pruning_strategy: str = "none", 
        pruning_parameter: str = None,          
        projection: str = "none", 
        predict_self = False,   
        time_strategy: str = "steady_state",  
        do_parallel: bool = True,
        kwargs = None,
    ):
        """Fit the model.

        Args:
            method (str): Regression method to use. Defaults to "RidgeCVExtraPenalty", which uses 
                sklearn.linear_model.RidgeCV and combats overfitting by scanning higher penalty params whenever
                the highest one is selected. Other methods not implemented yet. 
            confounders (list): Not implemented yet.
            cell_type_sharing_strategy (str, optional): Whether to fit one model across all training data ('identical') 
                or fit separate ones ('distinct'). Defaults to "distinct".
            cell_type_labels (str): Name of column in self.train.obs to use for cell type labels.
            network_prior (str, optional): How to incorporate user-provided network structure. 
                - "ignore": don't use it. 
                - "restrictive" (default): allow only user-specified regulators for each target.
                - maybe more options will be implemented. 
            pruning_strategy (str, optional) 
                - "prune_and_refit": keep the largest n coefficients across all models (not per model), where n is pruning_parameter.
                - "none": don't prune the model.
                - maybe more options will be implemented. 
            pruning_parameter (numeric, optional): e.g. lasso penalty or total number of nonzero coefficients. See "pruning_strategy" for details.
            projection (str, optional): Not implemented yet.
            predict_self (bool, optional): Should e.g. POU5F1 activity be used to predict POU5F1 expression? Defaults to False.
            time_strategy (str): 'steady_state' predicts each a gene from sample i using TF activity features derived
                from sample i. 'two_step' trains a model to gradually transform control samples into perturbed samples by
                first perturbing the targeted gene, then propagating the perturbation to other TFs, then propagating throughout the genome. 
                Under development circa 2022-Dec-01; see Eric's slides for a cleaner explanation.
            kwargs: Passed to DCDFG. See help(dcdfg_wrapper.DCDFGWrapper.train). 
        """
        # Save some training args to later ensure prediction is consistent with training.
        self.training_args["confounders"]                = confounders
        self.training_args["network_prior"]              = network_prior
        self.training_args["cell_type_sharing_strategy"] = cell_type_sharing_strategy
        self.training_args["cell_type_labels"]           = cell_type_labels
        self.training_args["predict_self"]               = predict_self
        self.training_args["method"]                     = method
        self.training_args["time_strategy"]              = time_strategy
        if network_prior != "ignore" and self.training_args["cell_type_sharing_strategy"] != 'identical':
            ct_from_network = self.network.get_all_one_field("cell_type")
            ct_from_trainset = self.train.obs[cell_type_labels]
            ct_missing = [ct for ct in ct_from_trainset if ct not in ct_from_network]
            prettyprint = lambda x: '\n'.join(x)
            if len(ct_from_network) > 0:
                if len(ct_missing)>0:
                    raise ValueError(
                        "Some cell types in the training data are not in the networks."
                        "Trainset: \n" f"{prettyprint(ct_from_trainset)}"
                        "Networks: \n" f"{prettyprint(ct_from_network)}"
                        "Missing:  \n" f"{prettyprint(ct_missing)}"
                    )
        if method.startswith("DCDFG"):
            assert len(confounders)==0, "DCDFG cannot currently include confounders."
            assert network_prior=="ignore", "DCDFG cannot currently include known network structure."
            assert cell_type_sharing_strategy=="identical", "DCDFG cannot currently fit each cell type separately."
            assert not predict_self, "DCDFG cannot include autoregulation."
            assert time_strategy == 'steady_state', "DCDFG assumes steady state."
            factor_graph_model = dcdfg_wrapper.DCDFGWrapper()
            _, constraint_mode, model_type, do_use_polynomials = method.split("-")
            print(f"""DCDFG args parsed as:
               constraint_mode: {constraint_mode}
                    model_type: {model_type}
            do_use_polynomials: {do_use_polynomials}
            """)
            do_use_polynomials = do_use_polynomials =="True"
            self.models = factor_graph_model.train(
                self.train,
                constraint_mode = constraint_mode,
                model_type = model_type,
                do_use_polynomials = do_use_polynomials,
                **kwargs,
            )
            return 
        elif method == "mean":
            def FUN(X,y):
                return sklearn.dummy.DummyRegressor(strategy="mean").fit(X, y)
        elif method == "median":
            def FUN(X,y):
                return sklearn.dummy.DummyRegressor(strategy="median").fit(X, y)
        elif method == "GradientBoostingRegressor":
            def FUN(X,y):
                return sklearn.ensemble.GradientBoostingRegressor().fit(X, y)
        elif method == "ExtraTreesRegressor":
            def FUN(X,y):
                return sklearn.ensemble.ExtraTreesRegressor().fit(X, y)
        elif method == "KernelRidge":
            def FUN(X,y):
                return sklearn.kernel_ridge.KernelRidge().fit(X, y)
        elif method == "ElasticNetCV":
            def FUN(X,y):
                return sklearn.linear_model.ElasticNetCV(
                    fit_intercept=True,
                ).fit(X, y)
        elif method == "LarsCV":
            def FUN(X,y):
                return sklearn.linear_model.LarsCV(
                    fit_intercept=True,
                    normalize=False,
                ).fit(X, y)
        elif method == "OrthogonalMatchingPursuitCV":
            def FUN(X,y):
                return sklearn.linear_model.OrthogonalMatchingPursuitCV(
                    fit_intercept=True, 
                    normalize=False,
                ).fit(X, y)
        elif method == "ARDRegression":
            def FUN(X,y):
                return sklearn.linear_model.ARDRegression(
                    fit_intercept=True,
                ).fit(X, y)
        elif method == "BayesianRidge":
            def FUN(X,y):
                return sklearn.linear_model.BayesianRidge(
                    fit_intercept=True,
                ).fit(X, y)
        elif method == "LassoCV":
            def FUN(X,y):
                return sklearn.linear_model.LassoCV(
                    fit_intercept=True,
                ).fit(X, y)
        elif method == "LassoLarsIC":
            def FUN(X,y):
                return sklearn.linear_model.LassoLarsIC(
                    fit_intercept=True,
                    normalize=False,
                ).fit(X, y)
        elif method == "RidgeCV":
            def FUN(X,y):
                return sklearn.linear_model.RidgeCV(
                    alphas=(0.01, 0.1, 1.0, 10.0, 100), 
                    fit_intercept=True,
                    alpha_per_target=False, 
                    store_cv_values=True, #this lets us use ._cv_values later for simulating data.
                ).fit(X, y)
        elif method == "RidgeCVExtraPenalty":
            def FUN(X,y):
                rcv = sklearn.linear_model.RidgeCV(
                    alphas=(0.01, 0.1, 1.0, 10.0, 100), 
                    fit_intercept=True,
                    alpha_per_target=False, 
                    store_cv_values=True, #this lets us use ._cv_values later for simulating data.
                ).fit(X, y)
                if rcv.alpha_ == np.max(rcv.alphas):
                    bigger_alphas = rcv.alpha_ * np.array([0.1, 1, 10, 100, 1000, 10000, 100000])
                    rcv = sklearn.linear_model.RidgeCV(
                        alphas=bigger_alphas, 
                        fit_intercept=True,
                        alpha_per_target=False, 
                        store_cv_values=True,
                    ).fit(X, y)
                return rcv
        else:
            raise NotImplementedError(f"Method {method} is not supported.")
        self.apply_supervised_ml(
            FUN, 
            network_prior=network_prior, 
            pruning_strategy=pruning_strategy, 
            pruning_parameter=pruning_parameter,
            do_parallel=do_parallel,
        )

    def predict(
        self,
        perturbations: list,
        starting_expression: anndata.AnnData = None,
        do_parallel: bool = True,
        add_noise = False,
        noise_sd = None,
        seed = 0,
    ):
        """Predict expression after new perturbations.

        Args:
            perturbations (iterable of tuples): Iterable of tuples with gene and its expression after 
                perturbation, e.g. {("POU5F1", 0.0), ("NANOG", 0.0), ("non_targeting", np.nan)}. Anything with
                expression np.nan will be treated as a control, no matter the name.
            held_out_data (anndata.AnnData): Initial conditions in the same shape as the output predictions. If 
                None, this will be set to the mean of the training data control expression values.
            do_parallel (bool): if True, use joblib parallelization. 
            add_noise (bool): if True, return simulated data Y + e instead of predictions Y 
                where e is IID Gaussian with variance equal to the estimated residual variance.
            noise_sd (bool): sd of the variable e described above. Defaults to estimates from the fitted models.
            seed (int): RNG seed.
        """
        if starting_expression is None:
            starting_states = self.train.obs["is_control"].copy()
            starting_features   = self.features[starting_states,:].mean(axis=0, keepdims = True)


        # Set up containers for output & features
        nrow = len(perturbations)
        features = np.zeros((nrow, len(self.tf_list)))
        columns_to_transfer = self.training_args["confounders"].copy()
        if self.training_args["cell_type_sharing_strategy"] != "identical":
            columns_to_transfer.append(self.training_args["cell_type_labels"])
        predictions = anndata.AnnData(
            X = np.zeros((nrow, len(self.train.var_names))),
            dtype=np.float32,
            var = self.train.var.copy(),
            obs = pd.DataFrame(
                {
                    "perturbation":"NA", 
                    "expression_level_after_perturbation": -999
                }, 
                index = [str(i) for i in range(nrow)],
                columns = ["perturbation", "expression_level_after_perturbation"] + columns_to_transfer,
            )
        )
        # implement perturbations
        for i in range(len(perturbations)):
            idx_str = str(i)
            features[i, :] = starting_features
            # Expected input: comma-separated strings like ("C6orf226,TIMM50,NANOG", "0,0,0") 
            pert_genes = perturbations[i][0].split(",")
            pert_exprs = [float(f) for f in str(perturbations[i][1]).split(",")]
            assert len(pert_genes) == len(pert_exprs), f"Malformed perturbation in sample {i}: {perturbations[i][0]}, {perturbations[i][1]}"
            for pert_idx in range(len(pert_genes)):
                # If it's nan, leave it unperturbed -- used for studying fitted values on controls. 
                if not np.isnan(pert_exprs[pert_idx]):
                    column = [tf == pert_genes[pert_idx] for tf in self.tf_list]
                    features[i, column] = pert_exprs[pert_idx]
            predictions.obs.loc[idx_str, "perturbation"]                        = perturbations[i][0]
            predictions.obs.loc[idx_str, "expression_level_after_perturbation"] = perturbations[i][1]
            # Later on, we can make matched_control_idx more flexible, e.g. to control for covariates.
            matched_control_idx = np.where(self.train.obs["is_control"])[0]
            for col in columns_to_transfer:
                predictions.obs.loc[idx_str, col] = self.train.obs.loc[matched_control_idx, col].values

        # Make predictions
        if self.training_args["method"].startswith("DCDFG"):
            predictions = self.models.predict(perturbations)     
        else:
            if self.training_args["cell_type_sharing_strategy"] == "distinct":
                cell_type_labels = predictions.obs[self.training_args["cell_type_labels"]]
            else:
                cell_type_labels = None
            if do_parallel:
                def do_loop(network = self.network, features=features):
                    return Parallel(n_jobs=cpu_count()-1, verbose = 1, backend="loky")(
                        delayed(predict_one_gene)(
                            target = self.train.var_names[i],
                            model = self.models[i],
                            features = features, 
                            network = network,
                            training_args = self.training_args,
                            tf_list = self.tf_list,
                            cell_type_labels = cell_type_labels,                     
                        )
                        for i in range(len(self.train.var_names))
                    )
            else:
                def do_loop(network = self.network, features=features):
                    return [
                        predict_one_gene(
                            target = self.train.var_names[i],
                            model = self.models[i],
                            features = features, 
                            network = self.network,
                            training_args = self.training_args,
                            tf_list = self.tf_list,
                            cell_type_labels = cell_type_labels, 
                        )
                        for i in range(len(self.train.var_names))
                    ]
            y = do_loop()
            for i in range(len(self.train.var_names)):
                predictions.X[:,i] = y[i]
            # Set perturbed genes equal to user-specified expression, not whatever the endogenous level is predicted to be
            for i, pp in enumerate(perturbations):
                if pp[0] in predictions.var_names:
                    predictions[i, pp[0]].X = pp[1]
            # Do one more time-step
            if self.training_args["time_strategy"] == "two_step":
                y = do_loop(
                    features = self.extract_tf_activity(
                        train = predictions, 
                        in_place = False,
                        method = "tf_rna",
                    )
                )
                for i in range(len(self.train.var_names)):
                    predictions.X[:,i] = y[i]
                # Set perturbed genes equal to user-specified expression, not whatever the endogenous level is predicted to be
                for i, pp in enumerate(perturbations):
                    if pp[0] in predictions.var_names:
                        predictions[i, pp[0]].X = pp[1]
            
        # Add noise. This is useful for simulations. 
        if add_noise:
            np.random.seed(seed)
            for i in range(len(self.train.var_names)):
                if noise_sd is None:
                    # From the docs on the RidgeCV attribute ".cv_values_":
                    # 
                    # > only available if store_cv_values=True and cv=None). 
                    # > After fit() has been called, this attribute will contain 
                    # > the mean squared errors if scoring is None otherwise it will 
                    # > contain standardized per point prediction values.
                    # 
                    try:
                        noise_sd = np.sqrt(np.mean(self.models[i].cv_values_))
                    except AttributeError:
                        raise ValueError("Noise standard deviation could not be extracted from trained models. Please provide it when calling GRN.simulate().")
                predictions.X[:,i] = predictions.X[:,i] + np.random.standard_normal(len(predictions.X[:,i]))*noise_sd

        return predictions

# Having stand-alone functions (rather than methods of the GRN class) speeds
# up parallel computing by a FUCK TON, probably because it avoids copying the 
# entire GRN object to each new thread.  

def predict_one_gene(
    target: str, 
    model,
    features, 
    network,
    training_args,
    tf_list,
    cell_type_labels,
) -> np.array:
    """Predict expression of one gene after perturbation. 

    Args:
        model: self.models[i], but for parallelization, we avoid passing self because it would have to get serialized and sent to another process.
        train_obs: self.train.obs, but we avoid passing self because it would have to get serialized & sent to another process.
        network: self.network (usually), but we avoid passing self because it would have to get serialized & sent to another process.
        training_args: self.training_args, but we avoid passing self because it would have to get serialized & sent to another process.
        tf_list: self.tf_list, but we avoid passing self because it would have to get serialized & sent to another process.

        target (str): Which gene to predict.
        features (matrix-like): inputs to predictive model, with shape num_regulators by num_perturbations. 
        cell_type_labels (Iterable): in case models are cell type-specific, this specifies the cell types. 

    Returns:
        np.array: Expression of gene i across all perturbed samples. 
    """
    predictions = np.zeros(features.shape[0])
    if training_args["cell_type_sharing_strategy"] == "distinct":
        for cell_type in cell_type_labels.unique():
            is_in_cell_type = cell_type_labels==cell_type
            regulators = get_regulators(
                tf_list = tf_list, 
                predict_self = training_args["predict_self"], 
                network_prior = training_args["network_prior"],
                target = target, 
                network = network,
                cell_type=cell_type
            )
            is_in_model = [tf in regulators for tf in tf_list]
            if not any(is_in_model):
                X = 1 + 0*features[:,[0]]
            else:
                X = features[:,is_in_model]
            X = X[is_in_cell_type,:]
            predictions[is_in_cell_type] = model[cell_type].predict(X = X).squeeze()
    elif training_args["cell_type_sharing_strategy"] == "identical":
        regulators = get_regulators(
            tf_list = tf_list, 
            predict_self = training_args["predict_self"], 
            network_prior = training_args["network_prior"],
            target = target, 
            network = network,
            cell_type=None
        )
        is_in_model = [tf in regulators for tf in tf_list]
        if not any(is_in_model):
            X = 1 + 0*features[:,[0]]
        else:
            X = features[:,is_in_model]
        predictions = model.predict(X = X).squeeze()
        
    else:
        raise NotImplementedError("Invalid cell_type_sharing_strategy.")
    return predictions

def apply_supervised_ml_one_gene(
    train_obs,
    features,
    network,
    training_args, 
    tf_list, 
    FUN, 
    network_prior: str,
    target: str,
    target_expr
    ):
    """Apply a supervised ML method to predict target expression from TF activity (one target gene only).

    Args:

        train_obs: self.train.obs, but we avoid passing self because it would have to get serialized & sent to another process.
        features: self.features, but we avoid passing self so that we can treat this as memory-mapped. 
        network: self.network (usually), but we avoid passing self because it would have to get serialized & sent to another process.
        training_args: self.training_args, but we avoid passing self because it would have to get serialized & sent to another process.
        tf_list: self.tf_list, but we avoid passing self because it would have to get serialized & sent to another process.

        FUN (Callable): See GRN.apply_supervised_ml docs.
        network_prior (str): see GRN.fit docs.
        target (str): target gene symbol
        target_expr: target gene expression levels

    Returns:
        _type_: dict with cell types as keys and with values containing result of FUN. see GRN.fit docs.
    """
    is_target_perturbed = train_obs["perturbation"]==target
    if training_args["cell_type_sharing_strategy"] == "distinct":
        assert training_args["cell_type_labels"] is not None,               "cell_type_labels must be provided."
        assert training_args["cell_type_labels"] in set(train_obs.columns), "cell_type_labels must name a column in .obs of training data."
        models = {}
        for cell_type in train_obs[training_args["cell_type_labels"]].unique():
            rrrelevant_rrregulators = get_regulators(
                tf_list       = tf_list, 
                predict_self  = training_args["predict_self"], 
                network_prior = network_prior,
                target        = target, 
                network       = network,
                cell_type     = cell_type, 
            )
            is_in_model = [tf in rrrelevant_rrregulators for tf in tf_list]
            is_in_cell_type = train_obs[training_args["cell_type_labels"]]==cell_type
            if not any(is_in_model):
                X = np.ones(shape=features[:,[0]].shape)
            else:
                X = features[:,is_in_model]
            models[cell_type] = FUN(
                X = X[is_in_cell_type & ~is_target_perturbed,:], 
                y = target_expr[is_in_cell_type & ~is_target_perturbed]
            )
        return models
    elif training_args["cell_type_sharing_strategy"] == "similar":
        raise NotImplementedError("cell_type_sharing_strategy 'similar' is not implemented yet.")
    elif training_args["cell_type_sharing_strategy"] == "identical":
        rrrelevant_rrregulators = get_regulators(
            tf_list       = tf_list, 
            predict_self  = training_args["predict_self"], 
            network_prior = network_prior,
            target        = target, 
            network       = network,
            cell_type = None,
        )
        is_in_model = [tf in rrrelevant_rrregulators for tf in tf_list]
        if not any(is_in_model):
            X = np.ones(shape=features[:,[0]].shape)
        else:
            X = features[:,is_in_model]
        return FUN(X = X[~is_target_perturbed, :], y = target_expr[~is_target_perturbed])
    else:
        raise ValueError("cell_type_sharing_strategy must be 'distinct' or 'identical' or 'similar'.")
    return

def get_regulators(tf_list, predict_self, network_prior: str, target: str, network = None, cell_type = None) -> list:
    """Get candidates for what's directly upstream of a given gene.

    Args:
        tf_list: list of all candidate regulators
        predict_self: Should a candidate regulator be used to predict itself? 
        network_prior (str): see GRN.fit docs
        target (str): target gene
        network (load_networks.LightNetwork): see GRN.fit docs
        cell_type: if network structure is cell-type-specific, which cell type to use when getting regulators.

    Returns:
        (list of str): candidate regulators currently included in the model for predicting 'target'.
        Given the same inputs, these will always be returned in the same order -- furthermore, in the 
        order that they appear in tf_list.
    """
    if network_prior == "ignore":
        selected_features = tf_list.copy()
    elif network_prior == "restrictive":
        assert network is not None, "For restrictive network priors, you must provide the network as a LightNetwork object."
        regulators = network.get_regulators(target=target)
        if cell_type is not None and "cell_type" in regulators.keys():
            regulators = regulators.loc[regulators["cell_type"]==cell_type,:]
        selected_features = [tf for tf in tf_list if tf in set(regulators["regulator"])]
    else:
        raise ValueError("network_prior must be one of 'ignore' and 'restrictive'. ")

    if not predict_self:
        try:
            selected_features = [tf for tf in selected_features if tf != target]
        except (KeyError, ValueError) as e:
            pass
    if len(selected_features)==0:
        return ["NO_REGULATORS"]
    else:
        return selected_features

