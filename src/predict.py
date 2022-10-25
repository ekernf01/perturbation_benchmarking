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
os.chdir(PROJECT_PATH + "benchmarking")
DEFAULT_TF_LIST = pd.read_csv("../accessory_data/humanTFs.csv")
DEFAULT_TF_LIST = [g for g in DEFAULT_TF_LIST.loc[DEFAULT_TF_LIST["Is TF?"]=="Yes","HGNC symbol"]]
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbations', 'load_perturbations'))) 
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

    def get_regulators(self, network_prior: str, g: str):
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
            assert self.network is not None, "For restrictive network priors, you must provide the network as a pandas dataframe."
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

    def _apply_supervised_ml_one_gene(self, FUN, network_prior: str, g: str):
        """Apply a supervised ML method to predict target expression from TF activity (one target gene only).

        Args:
            FUN (Callable): See GRN.apply_supervised_ml docs.
            network_prior (str): see GRN.fit docs.
            g (str): target gene

        Returns:
            _type_: same as result of FUN. see GRN.fit docs.
        """
        rrrelevant_rrregulators = self.get_regulators(network_prior, g)
        is_in_model = [tf in rrrelevant_rrregulators for tf in self.tf_list]

        if not any(is_in_model):
            X = 1 + 0*self.features[:,[0]]
        else:
            X = self.features[:,is_in_model]

        return FUN(X = X, y = self.train[:,g].X)
        
    def apply_supervised_ml(self, FUN, network_prior, pruning_strategy, pruning_parameter, verbose:int = 1):
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
                delayed(self._apply_supervised_ml_one_gene)(  FUN, network_prior, g )
                for g in self.train.var_names
            )
            # Prune
            current_network = pd.concat(
                [
                    pd.DataFrame(
                        {
                            "regulator": self.get_regulators(network_prior, self.train.var_names[i]), 
                            "target": self.train.var_names[i], 
                            "weight": self.models[i].coef_.squeeze(),
                        } 
                    )
                    for i in range(len(self.train.var_names)) 
                ], 
                axis = 0)
            current_network = current_network.query("regulator != 'NO_REGULATORS'")
            current_network["abs_weight"] = np.abs(current_network["weight"])
            current_network = current_network.nlargest(pruning_parameter, "abs_weight")
            self.network = current_network
            # Refit
            self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose)(
                delayed(self._apply_supervised_ml_one_gene)(  FUN, network_prior, g )
                for g in self.train.var_names
            )
            
        elif pruning_strategy == "none":
            # self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose)(
            #     delayed(self._apply_supervised_ml_one_gene)(  FUN, network_prior, g )
            #     for g in self.train.var_names
            # )
            # This can help with debugging -- better backtrace than parallel code.
            self.models = [
                        self._apply_supervised_ml_one_gene(  FUN, network_prior, g )
                        for g in self.train.var_names
                    ]
        else:
            raise NotImplementedError("pruning_strategy should be one of 'none', 'lasso', 'prune_and_refit' ")
                

    def fit(
        self,
        method: str, 
        confounders: list, 
        cell_type_labels: str = None,
        cell_type_sharing: str = "distinct",   
        network_prior: str = "restrictive",    
        pruning_strategy: str = "none", 
        pruning_parameter: str = None,          
        projection: str = "none",      
    ):
        """Fit the model.

        Args:
            method (str): Regression method to use. Defaults to "linear", which uses sklearn.linear_model.RidgeCV. Others not implemented yet. 
            confounders (list): Not implemented yet.
            cell_type_labels (str): Not implemented yet.
            cell_type_sharing (str, optional): Not implemented yet. Defaults to "distinct".
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
        self.network_prior = network_prior
        if method == "linear":
            def FUN(X,y):
                return lm.RidgeCV(
                    alphas=(0.01, 0.1, 1.0, 10.0, 100), 
                    fit_intercept=True,
                    alpha_per_target=False
                ).fit(X, y)
        else:
            raise NotImplementedError("Only 'linear' is supported so far.")
        self.apply_supervised_ml(FUN, network_prior, pruning_strategy, pruning_parameter)

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
            adata.obs["perturbation"]                        = gene
            adata.obs["expression_level_after_perturbation"] = expression_level_after_perturbation
            return adata

        if starting_states is None:
            starting_states = self.train.obs["is_control"]
            
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
            regulators = self.get_regulators(self.network_prior, g)
            is_in_model = [g in regulators for g in self.tf_list]
            if not any(is_in_model):
                X = 1 + 0*features[:,[0]]
            else:
                X = features[:,is_in_model]
            predictions.X[:, i] = self.models[i].predict(X = X).squeeze()
        return predictions
