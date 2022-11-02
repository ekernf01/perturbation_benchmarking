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
            self.features = self.train[:,self.tf_list].X
        else:
            raise NotImplementedError("Only 'tf_rna' feature extraction is so far available.")
        
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
            self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose, prefer="threads")(
                delayed(apply_supervised_ml_one_gene)(
                    train_obs = self.train.obs,
                    target_expr = self.train[:,g].X,
                    features = self.features,
                    network = self.network,
                    training_args = self.training_args,
                    tf_list = self.tf_list, 
                    predict_self = self.predict_self, 
                    FUN=FUN, 
                    network_prior=network_prior,
                    target = g,
                )
                for g in self.train.var_names
            )
            # Prune
            chunks = []
            for cell_type in self.train.obs[self.training_args["cell_type_labels"]].unique():
                for i in range(len(self.train.var_names)):
                    chunks.append(
                        pd.DataFrame(
                                {
                                    "regulator": list(get_regulators(
                                            tf_list = self.tf_list, 
                                            predict_self = self.predict_self, 
                                            network_prior = network_prior,
                                            g = self.train.var_names[i], 
                                            network = self.network,
                                         )), 
                                    "target": self.train.var_names[i], 
                                    "weight": self.models[i][cell_type].coef_.squeeze(),
                                    "cell_type": cell_type,
                                } 
                            )
                        )
            pruned_network = pd.concat(chunks)
            pruned_network = pruned_network.query("regulator != 'NO_REGULATORS'")
            pruned_network.loc[:, "abs_weight"] = np.abs(pruned_network["weight"])
            pruned_network = pruned_network.groupby("cell_type")
            pruned_network = pruned_network.apply(lambda grp: grp.nlargest(pruning_parameter, "abs_weight"))
            # Refit
            self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose, prefer="threads")(
                delayed(apply_supervised_ml_one_gene)(  
                    train_obs = self.train.obs,
                    target_expr = self.train[:,g].X,
                    features = self.features,
                    network = pruned_network,
                    training_args = self.training_args,
                    tf_list = self.tf_list, 
                    predict_self = self.predict_self, 
                    FUN=FUN, 
                    network_prior=network_prior,
                    target = g,
                )                
                for g in self.train.var_names
            )
            
        elif pruning_strategy == "none":
            self.models = Parallel(n_jobs=cpu_count()-1, verbose = verbose, prefer="threads")(
                delayed(apply_supervised_ml_one_gene)(  
                    train_obs = self.train.obs,
                    target_expr = self.train[:,g].X,
                    features = self.features,
                    network = self.network,
                    training_args = self.training_args,
                    tf_list = self.tf_list, 
                    predict_self = self.predict_self, 
                    FUN=FUN, 
                    network_prior=network_prior,
                    target = g,
                )                
                for g in self.train.var_names
            )
            # # This can help with debugging -- better backtrace than you get with parallel code.
            # self.models = [
            #   apply_supervised_ml_one_gene(    
                    # train = self.train,
                    # features = self.features,
                    # network = self.network,
                    # training_args = self.training_args,
                    # FUN=FUN, 
                    # network_prior=network_prior,
                    # target = g,            
            #      )
            #     for g in self.train.var_names
            # ]
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

        def perturb_features(features, gene, expression_level_after_perturbation):
            features[:, self.tf_list == gene] = expression_level_after_perturbation
            return features
        def perturb_metadata(adata, gene, expression_level_after_perturbation):
            adata.obs.loc[:,"perturbation"]                        = gene
            adata.obs.loc[:,"expression_level_after_perturbation"] = expression_level_after_perturbation
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
        for i in range(len(self.models)):
            g = self.train.var_names[i]
            for cell_type in      predictions.obs[self.training_args["cell_type_labels"]].unique():
                is_in_cell_type = predictions.obs[self.training_args["cell_type_labels"]]==cell_type
                regulators = get_regulators(
                    tf_list = self.tf_list, 
                    predict_self = self.predict_self, 
                    network_prior = self.training_args["network_prior"],
                    g = g, 
                    network = self.network,
                )
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
        assert training_args["cell_type_labels"] in train_obs.columns, "cell_type_labels must name a column in .obs of training data."
        assert pd.api.types.is_string_dtype(train_obs[training_args["cell_type_labels"]]), "cell_type_labels must name a string column in .obs of training data."
    elif training_args["cell_type_sharing_strategy"] == "similar":
        raise NotImplementedError("cell_type_sharing_strategy 'similar' is not implemented yet.")
    elif training_args["cell_type_sharing_strategy"] == "identical":
        assert training_args["cell_type_labels"] == "one_giant_cell_type"
        assert all(train_obs.loc[:,training_args["cell_type_labels"]] == "all")
    else:
        raise ValueError("cell_type_sharing_strategy must be 'distinct' or 'identical' or 'similar'.")
    
    models = {}
    for cell_type in train_obs[training_args["cell_type_labels"]].unique():
        if network is not None and "cell_type" in network.keys():
            cell_type_subnetwork = network.loc[network["cell_type"]==cell_type,:]
        else:
            cell_type_subnetwork = network
        rrrelevant_rrregulators = get_regulators(
                                            tf_list       = tf_list, 
                                            predict_self  = predict_self, 
                                            network_prior = network_prior,
                                            g             = target, 
                                            network       = cell_type_subnetwork,
            )
        is_in_model = [tf in rrrelevant_rrregulators for tf in tf_list]
        is_in_cell_type = train_obs[training_args["cell_type_labels"]]==cell_type
        if not any(is_in_model):
            X = np.ones(shape=features[:,[0]].shape)
        else:
            X = features[:,is_in_model]
        models[cell_type] = FUN(X = X[is_in_cell_type,:], y = target_expr[is_in_cell_type])
    return models

def get_regulators(tf_list, predict_self, network_prior: str, g: str, network = None):
    """_summary_

    Args:
        tf_list: list of all candidate regulators
        predict_self: Should a candidate regulator be used to predict itself? 
        network_prior (str): see GRN.fit docs
        g (str): target gene

    Returns:
        (list of str): candidate regulators currently included in the model for predicting g.
    """
    if network_prior == "ignore":
        selected_features = tf_list.copy()
    elif network_prior == "restrictive":
        assert network is not None, "For restrictive network priors, you must provide the network as a pandas dataframe."
        selected_features = set(network.loc[network["target"]==g,"regulator"]).intersection(tf_list)
    else:
        raise ValueError("network_prior must be one of 'ignore' and 'restrictive'. ")

    if not predict_self:
        try:
            selected_features.remove(g)
        except (KeyError, ValueError) as e:
            pass
    if len(selected_features)==0:
        return ["NO_REGULATORS"]
    else:
        return selected_features