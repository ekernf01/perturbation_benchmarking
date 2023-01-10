import torch
import pytorch_lightning as pl
import ggrn_backend3.linear
import anndata
import math
import pandas as pd
import numpy as np
from pytorch_lightning.callbacks.early_stopping import EarlyStopping

class AnnDataMatchedControlsDataSet(torch.utils.data.Dataset):
    def __init__(self, adata: anndata.AnnData, matching_method: str) -> None:
        super().__init__()
        self.adata = adata
        self.adata.var["numeric_index"] = [i for i in range(self.adata.var.shape[0])]
        def apply_over_comma_separated_string(s: str, f, sep=","):
            return sep.join([str(f(x)) for x in s.split(sep)])
        self.adata.obs["perturbation_index"] = [
            apply_over_comma_separated_string(p, self.get_index_from_gene)
            for p in self.adata.obs["perturbation"]
        ]

        # Mark each observation as treatment, control, and/or steady state
        if "is_control" not in set(self.adata.obs.columns):
            raise KeyError("train data must have a boolean column in .obs with key 'is_control'.")
        if self.adata.obs["is_control"].dtype.name != "bool":
            dt = self.adata.obs["is_control"].dtype
            raise TypeError(f"train.obs['is_control'] must be boolean. dtype: {repr(dt)}")
        if "is_treatment" not in set(self.adata.obs.columns):
            self.adata.obs["is_treatment"] = ~self.adata.obs["is_control"]
        if "is_steady_state" not in set(self.adata.obs.columns):
            self.adata.obs["is_steady_state"] = self.adata.obs["is_control"]
        if "perturbation" not in self.adata.obs.columns:
            raise KeyError("train data must have a comma-separated str column in .obs with key 'perturbation'.")
        if "expression_level_after_perturbation" not in self.adata.obs.columns:
            raise KeyError("train data must have a comma-separated str column in .obs with key 'expression_level_after_perturbation'.")
        if self.adata.obs["expression_level_after_perturbation"].dtype.name != "str":
            self.adata.obs['expression_level_after_perturbation'] = self.adata.obs['expression_level_after_perturbation'].astype(str)
        self.adata = MatchControls(self.adata, matching_method)

    def __len__(self):
        return self.adata.obs.shape[0]

    def __getitem__(self, idx):
        matched_idx = self.adata.obs["matched_control"][idx]
        return {
            "treatment":      {
                "expression": self.adata.X[idx,        :], 
                "metadata": self.adata.obs.iloc[idx,:].to_dict(),
            },
            "matched_control":{
                "expression":        self.adata.X[matched_idx,:],
                "metadata":   self.adata.obs.iloc[matched_idx,:].to_dict(),
            },
        }
    
    def get_index_from_gene(self, g):
        return self.adata.var.loc[g, "numeric_index"] if g in self.adata.var_names else -999 #sentinel value indicates this is a control

def MatchControls(train_data: anndata.AnnData, matching_method: str):
    if matching_method == "closest":
        raise NotImplementedError("Cannot yet match to closest control.")
        train_data.obs["matched_control"] = "placeholder"
    elif matching_method == "random":
        train_data.obs["matched_control"] = np.random.choice(
            np.where(train_data.obs["is_control"])[0], 
            train_data.obs.shape[0], 
            replace = True,
        )
    elif matching_method == "user":
        assert "matched_control" in train_data.obs.columns
    else: 
        raise ValueError("matching method must be 'closest' or 'random'.")
    return train_data

def kaiming_init(model):
    """
    Initialization suitable for use with leaky ReLU, from Kaiming He et al 2015.
    From the official Pytorch website, https://pytorch-lightning.readthedocs.io/en/stable/notebooks/course_UvA-DL/03-initialization-and-optimization.html
    Author: Phillip Lippe
    License: CC-by-SA
    """
    for name, param in model.named_parameters():
        if name.endswith(".bias"):
            param.data.fill_(0)
        elif name.startswith("layers.0"):  # The first layer does not have ReLU applied on its input
            param.data.normal_(0, 1 / math.sqrt(param.shape[1]))
        else:
            param.data.normal_(0, math.sqrt(2) / math.sqrt(param.shape[1]))

class GGRNAutoregressiveModel:

    def __init__(self, train_data, matching_method):
        self.train_data = AnnDataMatchedControlsDataSet(train_data, matching_method)
        self.model = None
        return

    def train(
        self,
        S = 1,
        regression_method = "linear",
        low_dimensional_structure = "none",
        low_dimensional_training = "PCA",
        network = None,
        device = None,
        limit_train_batches = 100,
        max_epochs = 1, 
        learning_rate = 0.001,
        batch_size = 64,
        regularization_parameter = 0,
    ):
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        self.model = ggrn_backend3.linear.LinearAutoregressive(
            n_genes = self.train_data.adata.X.shape[1],
            S = S,
            regression_method = regression_method,
            low_dimensional_structure = low_dimensional_structure,
            low_dimensional_training = low_dimensional_training,
            learning_rate = learning_rate, 
            regularization_parameter = regularization_parameter,
        )
        kaiming_init(self.model)

        dataloader = torch.utils.data.DataLoader(self.train_data, batch_size=batch_size, shuffle=True)
        trainer = pl.Trainer(
            limit_train_batches=limit_train_batches, 
            max_epochs=max_epochs,
            accelerator=device,
            callbacks=[EarlyStopping(monitor="training_loss", mode="min")]
        )

        trainer.fit(model=self.model, train_dataloaders=dataloader)
        return 
        
    def predict(
        self, 
        example_perturbations, 
        ):
        # TODO: select controls
        # TODO: run model forwards
        # self.model()
        return

