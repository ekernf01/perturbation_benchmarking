import argparse
import os
import anndata
import pandas as pd

import numpy as np
import pytorch_lightning as pl
from pytorch_lightning.callbacks import EarlyStopping
from pytorch_lightning.loggers import WandbLogger
from torch.utils.data import DataLoader, random_split

import wandb
import sys
os.chdir("/home/gary/cahan_rotation/perturbation_benchmarking/src/dcdfg_wrapper")
sys.path.append("/home/gary/cahan_rotation/perturbation_benchmarking/src/dcdfg_wrapper")
from dcdfg.callback import (AugLagrangianCallback, ConditionalEarlyStopping,
                            CustomProgressBar)
from dcdfg.linear_baseline.model import LinearGaussianModel
from dcdfg.lowrank_linear_baseline.model import LinearModuleGaussianModel
from dcdfg.lowrank_mlp.model import MLPModuleGaussianModel
from dcdfg.perturbseq_data import PerturbSeqDataset

import torch
torch.set_default_tensor_type(torch.DoubleTensor)

class DCDFGWrapper:
    def __init__(self):
        self.model = None
        self.train_dataset = None
    
    def train(
        self,
        adata: anndata.AnnData, 
        train_batch_size: int = 64,
        num_train_epochs: int = 600,
        num_fine_epochs: int = 50,
        num_modules: int = 20,
        learning_rate: float = 1e-3,
        regularization_parameter: float = 0.1,
        constraint_mode: str = "spectral_radius",
        model_type: str = "linearlr",    
        do_use_polynomials: bool = False,
        logfile: str = "logs",
        num_gpus = 0,
    ):
        # TODO: figure out how to feed this an anndata object instead of a filename
        train_dataset = PerturbSeqDataset(adata)
        print("Train dataset size", train_dataset.dim, len(train_dataset))
        nb_nodes = train_dataset.dim
        train_size = int(0.8 * len(train_dataset))
        val_size = len(train_dataset) - train_size
        train_dataset, val_dataset = random_split(train_dataset, [train_size, val_size])
        self.train_dataset = train_dataset
        
        if model_type == "linear":
            # create model
            self.model = LinearGaussianModel(
                nb_nodes,
                lr_init=learning_rate,
                reg_coeff=regularization_parameter,
                constraint_mode=constraint_mode,
                poly=do_use_polynomials,
            )
        elif model_type == "linearlr":
            self.model = LinearModuleGaussianModel(
                nb_nodes,
                num_modules,
                lr_init=learning_rate,
                reg_coeff=regularization_parameter,
                constraint_mode=constraint_mode,
            )            
        elif model_type == "mlplr":
            self.model = MLPModuleGaussianModel(
                nb_nodes,
                2,
                num_modules,
                16,
                lr_init=learning_rate,
                reg_coeff=regularization_parameter,
                constraint_mode=constraint_mode,
            )
        elif model_type == "already_trained":
            assert self.model is not None
        else:
            raise ValueError("model_type should be linear or linearlr or mlplr or already_trained.")

        # LOG CONFIG
        logger = WandbLogger(project=logfile, log_model=True)
        model_name = self.model.__class__.__name__
        if do_use_polynomials and model_name == "LinearGaussianModel":
            model_name += "_poly"
        logger.experiment.config.update(
            {"model_name": model_name, "module_name": self.model.module.__class__.__name__}
        )
        
        # Step 1: augmented lagrangian
        early_stop_1_callback = ConditionalEarlyStopping(
            monitor="Val/aug_lagrangian",
            min_delta=1e-4,
            patience=5,
            verbose=True,
            mode="min",
        )
        trainer = pl.Trainer(
            gpus=num_gpus,
            max_epochs=num_train_epochs,
            logger=logger,
            val_check_interval=1.0,
            callbacks=[AugLagrangianCallback(), early_stop_1_callback, CustomProgressBar()],
        )
        trainer.fit(
            self.model,
            DataLoader(train_dataset, batch_size=train_batch_size, num_workers=4),
            DataLoader(val_dataset, num_workers=8, batch_size=256),
        )
        wandb.log({"nll_val": self.model.nlls_val[-1]})
        wandb.finish()

        # freeze and prune adjacency
        self.model.module.threshold()
        # WE NEED THIS BECAUSE IF it's exactly a DAG THE POWER ITERATIONS DOESN'T CONVERGE
        # TODO Just refactor and remove constraint at validation time
        self.model.module.constraint_mode = "exp"
        # remove dag constraints: we have a prediction problem now!
        self.model.gamma = 0.0
        self.model.mu = 0.0

        # Step 2:fine tune weights with frozen model
        logger = WandbLogger(project=logfile, log_model=True)
        model_name = self.model.__class__.__name__
        if do_use_polynomials and model_name == "LinearGaussianModel":
            model_name += "_poly"
        logger.experiment.config.update(
            {"model_name": model_name, "module_name": self.model.module.__class__.__name__}
        )

        early_stop_2_callback = EarlyStopping(
            monitor="Val/nll", min_delta=1e-6, patience=5, verbose=True, mode="min"
        )
        trainer_fine = pl.Trainer(
            gpus=num_gpus,
            max_epochs=num_fine_epochs,
            logger=logger,
            val_check_interval=1.0,
            callbacks=[early_stop_2_callback, CustomProgressBar()],
        )
        trainer_fine.fit(
            self.model,
            DataLoader(train_dataset, batch_size=train_batch_size),
            DataLoader(val_dataset, num_workers=2, batch_size=256),
        )

        # Original implementation extracts the network structure for later inspection
        pred_adj = self.model.module.weight_mask.detach().cpu().numpy()
        
        assert np.equal(np.mod(pred_adj, 1), 0).all()

        # Step 4: add valid nll and dump metrics
        pred = trainer_fine.predict(
            ckpt_path="best",
            dataloaders=DataLoader(val_dataset, num_workers=8, batch_size=256),
        )
        val_nll = np.mean([x.item() for x in pred])
    
        acyclic = int(self.model.module.check_acyclicity())
        if acyclic != 1:
            print("""
            Graph is not acyclic! Problems:
            
            - Joint density may not be well-defined, so NLL may not be a fair evaluation metric
            - Predicted expression may diverge as effects propagate in cycles (yields NaN's)
            
            Possible solutions:

            - Retrain with different inputs (e.g. more epochs)
            - Continue training where you left off by calling the 'train' method with model_type == "already_trained"
            - Avoid using likelihood, and avoid divergence by using a small 'maxiter' in the 'simulateKO()' method.
            
            """)
        wandb.log(
            {
                "val nll": val_nll,
                "acyclic": acyclic,
                "n_edges": pred_adj.sum(),
            }
        )
        
        wandb.finish()
        
        return self

    def predict(self, perturbations, baseline_expression = None):
        if baseline_expression is None:
            baseline_expression = self.train_dataset.dataset.adata.X[
                self.train_dataset.dataset.adata.obs["is_control"],:
            ].mean(axis=0)
        genes = self.train_dataset.dataset.adata.var_names
        predicted_adata = anndata.AnnData(
            X = np.zeros((len(perturbations), len(genes))),
            var = self.train_dataset.dataset.adata.var.copy(),
            obs = pd.DataFrame(
                {
                    "perturbation"                       :[p[0] for p in perturbations], 
                    "expression_level_after_perturbation":[p[1] for p in perturbations],
                },
                index = range(len(perturbations)), 
            ),
            dtype = np.float64
        )
        def convert_gene_symbol_to_index(gene):
            return [i for i,g in enumerate(genes) if g==gene]
        
        with torch.no_grad():
            predicted_adata.X = self.model.simulateKO(
                control_expression = baseline_expression,
                KO_gene_indices    = [convert_gene_symbol_to_index(pp[0]) for pp in perturbations],
                KO_gene_values     = [pp[1] for pp in perturbations],
                maxiter            = 100
            )
        return predicted_adata
