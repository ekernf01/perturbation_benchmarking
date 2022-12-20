import torch
import pytorch_lightning as pl
import linear

class AutoregressiveModel:
    def __init__():
        self.model = None
        return
    def train(
        train_data,
        is_treatment = None, # defaults to ~train_data.obs["is_control"]
        is_control   = None, # defaults to train_data.obs["is_control"]
        is_steady_state = None, # defaults to is_control
        matching = "closest",
        S = 1,
        regression_method = "linear",
        low_dimensional_structure = "none",
        low_dimensional_training = "PCA",
        network = None,
        device = None,
    ):
        if device is None:
            device = "cuda" if torch.cuda.is_available() else "cpu"
        # TODO: make dataloader
        # TODO: select appropriate model architecture
        self.model = linear.LinearAutoregressive()
        # TODO: run optimizer        
        return 
        
    def predict(example_perturbations):
        # TODO: select controls
        # TODO: run model forwards
        self.model()

