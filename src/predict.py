from multiprocessing.sharedctypes import Value
import anndata
from joblib import Parallel, delayed, cpu_count
import pandas as pd
import os
import sklearn.linear_model as lm
import numpy as np

# Project-specific paths
import sys 
import importlib
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
DEFAULT_TF_LIST = pd.read_csv("../accessory_data/humanTFs.csv")
DEFAULT_TF_LIST = [g for g in DEFAULT_TF_LIST.loc[DEFAULT_TF_LIST["Is TF?"]=="Yes","HGNC symbol"]]
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 
import load_perturbations
importlib.reload(load_perturbations)

class GRN:
    """
    Flexible inference of gene regulatory network models.
    """
    def __init__(self, train: anndata.AnnData, network: pd.DataFrame = None, tf_list = DEFAULT_TF_LIST, predict_self = False, validate_immediately = True):
        """Create a GRN object.

        Args:
            train (anndata.AnnData): Training data. Should conform to the requirements in load_perturbations.check_perturbation_dataset().
            network (pd.DataFrame, optional): Dataframe containing networks. This arg is experimental. 
                The format and meaning may change drastically. Currently this expects columns ['tf', 'target', 'weight', 'tag'].
            tf_list (_type_, optional): List of gene names that are allowed to be regulators.
            predict_self (bool, optional): Should e.g. POU5F1 activity be used to predict POU5F1 expression? Defaults to False.
        """
        self.train = train 
        self.network = network 
        self.tf_list = list(set(tf_list).intersection(set(train.var_names)))
        self.predict_self = predict_self
        self.models = []
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
            self.features = self.train[:,self.tf_list].X.copy()
        else:
            raise NotImplementedError("Only 'tf_rna' feature extraction is so far available.")

    def get_regulators(self, network_prior: str, g: str, network_structure = None):
        """_summary_

        Args:
            network_prior (str): see GRN.fit docs
            g (str): target gene

        Returns:
            (list of str): candidate regulators currently included in the model for predicting g.
        """
        assert g in self.train.var_names, "This target is not in the training data!"
        if network_prior == "ignore":
            selected_features = self.tf_list.copy()
        elif network_prior == "restrictive":
            if network_structure is None:
                network_structure = self.network
            assert network_structure is not None, "For restrictive network priors, you must provide the network as a pandas dataframe."
            selected_features = set(self.network.loc[self.network["target"]==g,"regulator"]).intersection(self.tf_list)
        else:
            raise ValueError("network_prior must be one of 'ignore' and 'restrictive'. ")

        if not self.predict_self:
            try:
                selected_features.remove(g)
            except (KeyError, ValueError) as e:
                pass
        if len(selected_features)==0:
            return ["NO_REGULATORS"]
        else:
            return selected_features

    def _apply_supervised_ml_one_gene(
        self, 
        FUN, 
        network_prior: str,
        target: str,
        network_structure = None,
        ):
        """Apply a supervised ML method to predict target expression from TF activity (one target gene only).

        Args:
            FUN (Callable): See GRN.apply_supervised_ml docs.
            network_prior (str): see GRN.fit docs.
            target (str): target gene

        Returns:
            _type_: dict with cell types as keys and with values containing result of FUN. see GRN.fit docs.
        """
        if network_structure is None:
            network_structure = self.network
        if self.training_args["cell_type_sharing_strategy"] == "distinct":            
            assert self.training_args["cell_type_labels"] in self.train.obs.columns, "cell_type_labels must name a column in .obs of training data."
            assert pd.api.types.is_string_dtype(self.train.obs[self.training_args["cell_type_labels"]]), "cell_type_labels must name a string column in .obs of training data."
        elif self.training_args["cell_type_sharing_strategy"] == "similar":
            raise NotImplementedError("cell_type_sharing_strategy 'similar' is not implemented yet.")
        elif self.training_args["cell_type_sharing_strategy"] == "identical":
            assert self.training_args["cell_type_labels"] == "one_giant_cell_type"
            assert all(self.train.obs.loc[:,self.training_args["cell_type_labels"]] == "all")
        else:
            raise ValueError("cell_type_sharing_strategy must be 'distinct' or 'identical' or 'similar'.")
        
        models = {}
        for cell_type in self.train.obs[self.training_args["cell_type_labels"]].unique():
            if "cell_type" in network_structure.keys():
                cell_type_subnetwork = network_structure.loc[network_structure["cell_type"]==cell_type,:]
            else:
                cell_type_subnetwork = network_structure
            rrrelevant_rrregulators = self.get_regulators(network_prior, target, cell_type_subnetwork)
            is_in_model = [tf in rrrelevant_rrregulators for tf in self.tf_list]
            is_in_cell_type = self.train.obs[self.training_args["cell_type_labels"]]==cell_type
            if not any(is_in_model):
                X = 1 + 0*self.features[:,[0]]
            else:
                X = self.features[:,is_in_model]
            models[cell_type] = FUN(X = X[is_in_cell_type,:], y = self.train[is_in_cell_type,target].X)
        return models
        
    def apply_supervised_ml(
        self, 
        FUN, 
        network_prior, 
        pruning_strategy, 
        pruning_parameter,
        verbose:int = 1
    ):
        """Apply a supervised ML method to predict target expression from TF activity.

        Args:
            FUN (function): should accept inputs X,y and return a fitted model, sklearn-style, with a .predict method.
            verbose (int): passes to joblib.Parallel. https://joblib.readthedocs.io/en/latest/parallel.html
        """
        if pruning_strategy == "lasso":
            raise NotImplementedError("lasso pruning not implemented yet.")
        elif pruning_strategy == "prune_and_refit":
            # Fit
            self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose)(
                delayed(self._apply_supervised_ml_one_gene)(  
                    FUN=FUN, 
                    network_prior=network_prior,
                    target = g,
                )
                for g in self.train.var_names
            )
            # Prune
            pruned_network = pd.concat(
                [
                    pd.concat(
                        [
                            pd.DataFrame(
                                {
                                    "regulator": self.get_regulators(network_prior, self.train.var_names[i]), 
                                    "target": self.train.var_names[i], 
                                    "weight": self.models[i][cell_type].coef_.squeeze(),
                                    "cell_type": cell_type,
                                } 
                            )
                            for i in range(len(self.train.var_names)) 
                        ], 
                        axis = 0
                    ) 
                    for cell_type in self.train.obs[self.training_args["cell_type_labels"]]
                ])
            pruned_network = pruned_network.query("regulator != 'NO_REGULATORS'")
            pruned_network.loc[:, "abs_weight"] = np.abs(pruned_network["weight"])
            pruned_network = pruned_network.groupby("cell_type")
            pruned_network = pruned_network.apply(lambda grp: grp.nlargest(pruning_parameter, "abs_weight"))
            # Refit
            self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose)(
                delayed(self._apply_supervised_ml_one_gene)(  
                    FUN=FUN, 
                    network_prior=network_prior,
                    target = g,
                    network_structure = pruned_network,
                )                
                for g in self.train.var_names
            )
            
        elif pruning_strategy == "none":
            # self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose)(
            #     delayed(self._apply_supervised_ml_one_gene)(  
            #         FUN=FUN, 
            #         network_prior=network_prior,
            #         target = g,
            #     )                
            #     for g in self.train.var_names
            # )
            # This can help with debugging -- better backtrace than you get with parallel code.
            self.models = [
                self._apply_supervised_ml_one_gene(    
                    FUN=FUN, 
                    network_prior=network_prior,
                    target = g,
                 )
                for g in self.train.var_names
            ]
        else:
            raise NotImplementedError("pruning_strategy should be one of 'none', 'lasso', 'prune_and_refit' ")
                

    def fit(
        self,
        method: str, 
        confounders: list, 
        cell_type_sharing_strategy: str = "distinct",   
        cell_type_labels: str = None,
        network_prior: str = "restrictive",    
        pruning_strategy: str = "none", 
        pruning_parameter: str = None,          
        projection: str = "none",      
    ):
        """Fit the model.

        Args:
            method (str): Regression method to use. Defaults to "linear", which uses sklearn.linear_model.RidgeCV. Others not implemented yet. 
            confounders (list): Not implemented yet.
            cell_type_sharing_strategy (str, optional): Not implemented yet. Defaults to "distinct".
            cell_type_labels (str): Not implemented yet.
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
        self.training_args["network_prior"]              = network_prior
        self.training_args["cell_type_sharing_strategy"] = cell_type_sharing_strategy
        self.training_args["cell_type_labels"]           = cell_type_labels

        if self.training_args["cell_type_sharing_strategy"] == "identical":
            self.training_args["cell_type_labels"] = "one_giant_cell_type"
            self.train.obs.loc[:,self.training_args["cell_type_labels"]] = "all"

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
        )

    def predict(
        self,
        perturbations,
        starting_states = None,
    ):
        """Predict expression after new perturbations.

        Args:
            perturbations (iterable of tuples): Iterable of tuples with gene and its expression after 
                perturbation, e.g. {("POU5F1", 0.0), ("NANOG", 0.0)}.
            starting_states: indices of observations in self.train to use as initial conditions. Defaults to self.train.obs["is_control"].
        """
        def perturb_features(features, gene, expression_level_after_perturbation):
            features[:, self.tf_list == gene] = expression_level_after_perturbation
            return features
        def perturb_metadata(adata, gene, expression_level_after_perturbation):
            adata.obs.loc[:,"perturbation"]                        = gene
            adata.obs.loc[:,"expression_level_after_perturbation"] = expression_level_after_perturbation
            return adata

        if starting_states is None:
            starting_states = self.train.obs["is_control"].copy()
            
        features = np.concatenate([
            perturb_features(self.features[starting_states,:].copy(), p[0], p[1])
            for p in perturbations
        ])
        predictions = anndata.concat(
            perturb_metadata(self.train[   starting_states,:].copy(), p[0], p[1])
            for p in perturbations
        )
        for i in range(len(self.models)):
            g = self.train.var_names[i]
            for cell_type in      predictions.obs[self.training_args["cell_type_labels"]].unique():
                is_in_cell_type = predictions.obs[self.training_args["cell_type_labels"]]==cell_type
                regulators = self.get_regulators(self.training_args["network_prior"], g)
                is_in_model = [g in regulators for g in self.tf_list]
                if not any(is_in_model):
                    X = 1 + 0*features[:,[0]]
                else:
                    X = features[:,is_in_model]
                X = X[is_in_cell_type,:]
                M = self.models[i][cell_type]
                Y = M.predict(X = X).squeeze()
                predictions.X[is_in_cell_type, i] = Y
        return predictions
