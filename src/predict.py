from multiprocessing.sharedctypes import Value
import anndata
from joblib import Parallel, delayed, cpu_count
import pandas as pd
import os
import sklearn.linear_model as lm
import numpy as np
import gc 

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

class GRN:
    """
    Flexible inference of gene regulatory network models.
    """
    
    def __init__(
        self, train: anndata.AnnData, 
        network: load_networks.LightNetwork = None, 
        tf_list = DEFAULT_TF_LIST, 
        predict_self = False, 
        validate_immediately = True
        ):
        """Create a GRN object.

        Args:
            train (anndata.AnnData): Training data. Should conform to the requirements in load_perturbations.check_perturbation_dataset().
            network (pd.DataFrame, optional): LightNetwork object containing prior knowledge about regulators and targets.
            tf_list (_type_, optional): List of gene names that are allowed to be regulators.
            predict_self (bool, optional): Should e.g. POU5F1 activity be used to predict POU5F1 expression? Defaults to False.
        """
        self.train = train 
        assert network is None or type(network)==load_networks.LightNetwork
        self.network = network 
        self.tf_list = list(set(tf_list).intersection(set(train.var_names)))
        self.predict_self = predict_self
        self.models = [None for _ in self.train.var_names]
        self.training_args = {}
        if validate_immediately:
            assert self.check_perturbation_dataset()

    def check_perturbation_dataset(self):
        return load_perturbations.check_perturbation_dataset(ad=self.train)

    def extract_features(self, method = "tf_rna"):
        """Create a feature matrix where each row matches a row in self.train 
        each column represents activity of the corresponding TF in self.tf_list.

        Args:
            method (str, optional): How to extract features. Defaults to "tf_rna", which uses the mRNA level of the TF.
        """
        if method == "tf_rna":
            self.features = self.train[:,self.tf_list].X
        else:
            raise NotImplementedError("Only 'tf_rna' feature extraction is so far available.")

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
            _type_: _description_
        """
        if do_parallel:     
            def do_loop(network = self.network):
                    return Parallel(n_jobs=cpu_count()-1, verbose = verbose, backend="loky")(
                        delayed(apply_supervised_ml_one_gene)(
                            train_obs = self.train.obs,
                            target_expr = self.train[:,g].X,
                            features = self.features,
                            network = network,
                            training_args = self.training_args,
                            tf_list = self.tf_list, 
                            predict_self = self.predict_self, 
                            FUN=FUN, 
                            network_prior=network_prior,
                            target = g,
                        )
                        for g in self.train.var_names
                    )
        else:
            def do_loop(network = self.network):
                return [
                    apply_supervised_ml_one_gene(
                        train_obs = self.train.obs,
                        target_expr = self.train[:,g].X,
                        features = self.features,
                        network = network,
                        training_args = self.training_args,
                        tf_list = self.tf_list, 
                        predict_self = self.predict_self, 
                        FUN=FUN, 
                        network_prior=network_prior,
                        target = g,
                    )
                    for g in self.train.var_names
                ]

        if pruning_strategy == "lasso":
            raise NotImplementedError("lasso pruning not implemented yet.")
        elif pruning_strategy == "prune_and_refit":
            # Fit
            self.models = do_loop()
            # Prune
            chunks = []
            for i in range(len(self.train.var_names)):
                if self.training_args["cell_type_sharing_strategy"] == "distinct":
                    for cell_type in self.train.obs[self.training_args["cell_type_labels"]].unique():
                        chunks.append(
                            pd.DataFrame(
                                    {
                                        "regulator": get_regulators(
                                                tf_list = self.tf_list, 
                                                predict_self = self.predict_self, 
                                                network_prior = network_prior,
                                                target = self.train.var_names[i], 
                                                network = self.network,
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
                                            predict_self = self.predict_self, 
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
            # Refit
            self.models = do_loop(network = load_networks.LightNetwork(df=pruned_network))                
        elif pruning_strategy == "none":
            self.models = do_loop() 
        else:
            raise NotImplementedError("pruning_strategy should be one of 'none', 'lasso', 'prune_and_refit' ")
                
    
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
        do_parallel: bool = True,
    ):
        """Fit the model.

        Args:
            method (str): Regression method to use. Defaults to "linear", which uses sklearn.linear_model.RidgeCV. Others not implemented yet. 
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
        """
        # Save some training args to later ensure prediction is consistent with training.
        self.training_args["network_prior"]              = network_prior
        self.training_args["cell_type_sharing_strategy"] = cell_type_sharing_strategy
        self.training_args["cell_type_labels"]           = cell_type_labels

        if method == "linear":
            def FUN(X,y):
                return lm.RidgeCV(
                    alphas=(0.01, 0.1, 1.0, 10.0, 100), 
                    fit_intercept=True,
                    alpha_per_target=False
                ).fit(X, y)
        else:
            raise NotImplementedError("Only 'linear' is supported so far.")
        self.apply_supervised_ml(
            FUN, 
            network_prior=network_prior, 
            pruning_strategy=pruning_strategy, 
            pruning_parameter=pruning_parameter,
            do_parallel=do_parallel,
        )

    def predict(
        self,
        perturbations,
        starting_states = None,
        aggregation = "before",
    ):
        """Predict expression after new perturbations.

        Args:
            perturbations (iterable of tuples): Iterable of tuples with gene and its expression after 
                perturbation, e.g. {("POU5F1", 0.0), ("NANOG", 0.0)}.
            starting_states: indices of observations in self.train to use as initial conditions. Defaults to self.train.obs["is_control"].
            aggregation: one of "before", "after", or "none". If there are multiple starting states
        """

        if starting_states is None:
            starting_states = self.train.obs["is_control"].copy()
        
        # Handle aggregation of multiple controls
        if aggregation == "before":
            features          = self.features[starting_states,:].mean(axis=0, keepdims = True)
            metadata_template = self.train[starting_states,:][0,:]
        elif aggregation == "none":
            features = self.features[starting_states,:]
            metadata_template = self.train[starting_states,:]
        elif aggregation == "after":
            raise NotImplementedError("aggregation after simulation is not implemented yet, sorry.")
        else:
            raise ValueError("aggregation must be 'none', 'before', or 'after'.")

        # Initialize objects
        def perturb_features(features, gene, expression_level_after_perturbation):
            features[:, self.tf_list == gene] = expression_level_after_perturbation
            return features
        def perturb_metadata(adata, gene, expression_level_after_perturbation):
            adata.obs["perturbation"]                        = gene
            adata.obs["expression_level_after_perturbation"] = expression_level_after_perturbation
            return adata
        features = np.concatenate([
            perturb_features(features.copy(), p[0], p[1])
            for p in perturbations
        ])
        predictions = anndata.concat(
            keys = [
                p[0] + "_" + str(p[1])
                for p in perturbations
            ],
            adatas = [
                perturb_metadata(metadata_template.copy(), p[0], p[1])
                for p in perturbations
            ], 
            index_unique = "_",
        )
        # Make predictions
        for i in range(len(self.models)):
            target = self.train.var_names[i]
            if self.training_args["cell_type_sharing_strategy"] == "discrete":
                for cell_type in      predictions.obs[self.training_args["cell_type_labels"]].unique():
                    is_in_cell_type = predictions.obs[self.training_args["cell_type_labels"]]==cell_type
                    regulators = get_regulators(
                        tf_list = self.tf_list, 
                        predict_self = self.predict_self, 
                        network_prior = self.training_args["network_prior"],
                        target = target, 
                        network = self.network,
                        cell_type=cell_type
                    )
                    is_in_model = [tf in regulators for tf in self.tf_list]
                    if not any(is_in_model):
                        X = 1 + 0*features[:,[0]]
                    else:
                        X = features[:,is_in_model]
                    X = X[is_in_cell_type,:]
                    M = self.models[i][cell_type]
                    Y = M.predict(X = X).squeeze()
                    predictions.X[is_in_cell_type, i] = Y
            elif self.training_args["cell_type_sharing_strategy"] == "identical":
                regulators = get_regulators(
                    tf_list = self.tf_list, 
                    predict_self = self.predict_self, 
                    network_prior = self.training_args["network_prior"],
                    target = target, 
                    network = self.network,
                    cell_type=None
                )
                is_in_model = [tf in regulators for tf in self.tf_list]
                if not any(is_in_model):
                    X = 1 + 0*features[:,[0]]
                else:
                    X = features[:,is_in_model]
                M = self.models[i]
                Y = M.predict(X = X).squeeze()
                predictions.X[:, i] = Y
            else:
                raise NotImplementedError("Invalid cell_type_sharing_strategy.")
        return predictions

def apply_supervised_ml_one_gene(
    train_obs,
    features,
    network,
    training_args, 
    tf_list, 
    predict_self, 
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
        FUN (Callable): See GRN.apply_supervised_ml docs.
        network_prior (str): see GRN.fit docs.
        target (str): target gene symbol
        target_expr: target gene expression levels

    Returns:
        _type_: dict with cell types as keys and with values containing result of FUN. see GRN.fit docs.
    """
    if training_args["cell_type_sharing_strategy"] == "distinct":
        assert training_args["cell_type_labels"] is not None,          "cell_type_labels must name a column in .obs of training data."
        assert training_args["cell_type_labels"] in train_obs.columns, "cell_type_labels must name a column in .obs of training data."
        assert pd.api.types.is_string_dtype(train_obs[training_args["cell_type_labels"]]), "cell_type_labels must name a string column in .obs of training data."
        models = {}
        for cell_type in train_obs[training_args["cell_type_labels"]].unique():
            rrrelevant_rrregulators = get_regulators(
                tf_list       = tf_list, 
                predict_self  = predict_self, 
                network_prior = network_prior,
                target        = target, 
                network       = network,
                cell_type = cell_type,
            )
            is_in_model = [tf in rrrelevant_rrregulators for tf in tf_list]
            is_in_cell_type = train_obs[training_args["cell_type_labels"]]==cell_type
            if not any(is_in_model):
                X = np.ones(shape=features[:,[0]].shape)
            else:
                X = features[:,is_in_model]
            models[cell_type] = FUN(X = X[is_in_cell_type,:], y = target_expr[is_in_cell_type])
        return models
    elif training_args["cell_type_sharing_strategy"] == "similar":
        raise NotImplementedError("cell_type_sharing_strategy 'similar' is not implemented yet.")
    elif training_args["cell_type_sharing_strategy"] == "identical":
        rrrelevant_rrregulators = get_regulators(
            tf_list       = tf_list, 
            predict_self  = predict_self, 
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
        return FUN(X = X, y = target_expr)
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
        if cell_type is not None:
            assert "cell_type" in regulators.keys(), "'cell_type' must be a column in the provided network."
            regulators = regulators.loc[regulators["cell_type"]==cell_type,:]
        selected_features = [tf for tf in tf_list if tf in regulators["regulator"]]
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