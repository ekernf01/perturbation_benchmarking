import torch
import pytorch_lightning as pl
def P(x: torch.tensor, perturbations: list) -> torch.tensor:
    """Enact a perturbation.
    Args:
        x (torch.tensor): Gene expression vector
        p (list): list of tuples (target gene index, level), e.g. (5, 0) where gene 5 is Nanog and you knocked out Nanog.

    Returns:
        torch.tensor: Gene expression with targeted gene set to prescribed level. 
    """
    for p in perturbations:
        x[p[0]] = p[1]
    return x


class LinearAutoregressive(pl.LightningModule):
    def __init__(
        self, 
        n_genes,
        S,
        regression_method,
        low_dimensional_structure,
        low_dimensional_training,
        learning_rate, 
        regularization_parameter,
    ):
        super().__init__()
        self.S = S
        self.learning_rate = learning_rate
        self.regularization_parameter = regularization_parameter
        if regression_method == "linear":
            self.G = torch.nn.Linear(n_genes, n_genes)
        elif regression_method in {'multilayer_perceptron', 'mlp'}:
            raise NotImplementedError("mlp is not implemented yet.")
        else:
            raise ValueError(f"regression_method must be 'linear' or 'mlp'. Value: {regression_method}")
        
        if low_dimensional_structure in "none":
            if low_dimensional_training is not None:
                print("low_dimensional_training will be ignored since low_dimensional_structure is 'none'.")
        elif low_dimensional_structure in "RGQ":
            raise NotImplementedError("Low-D structure is not implemented yet.")
        else:
            raise ValueError("regression_method must be 'none' or 'RGQ'")

        if low_dimensional_training is None:
            pass
        elif low_dimensional_training == "supervised":
            raise NotImplementedError()
        elif low_dimensional_training == "PCA":
            raise NotImplementedError()
        elif low_dimensional_training == "fixed":
            raise NotImplementedError()
        else:
            raise ValueError(f"low_dimensional_training must be 'supervised','PCA', 'fixed'. Value: {low_dimensional_training}")

    def forward(self, x, perturbations):
        return self.G(P(x, perturbations))

    def training_step(self, input_batch):
        loss = 0
        for i in range(len(input_batch["treatment"]["metadata"]["perturbation"])):
            perturbations = zip(
                [int(g)   for g in input_batch["treatment"]["metadata"]["perturbation_index"][i].split(",")],
                [float(x) for x in input_batch["treatment"]["metadata"]["expression_level_after_perturbation"][i].split(",")],
            )
            perturbations = [p for p in perturbations if p[0]!=-999] #sentinel value indicates this is a control
            if input_batch["treatment"]["metadata"]["is_steady_state"][i]: 
                # (electric slide voice) one hop this time. *BEWMP*
                x_t = input_batch["treatment"]["expression"][i].clone()
                x_t = self(x_t, perturbations)
                loss += torch.linalg.norm(input_batch["treatment"]["expression"][i] - x_t)
            if input_batch["treatment"]["metadata"]["is_treatment"][i]: 
                # (electric slide voice) S hops this time. *BEWMP* *BEWMP* *BEWMP* *BEWMP* 
                x_t = input_batch["matched_control"]["expression"][i].clone()
                for _ in range(self.S):
                    x_t = self(x_t, perturbations)
                loss += torch.linalg.norm(input_batch["treatment"]["expression"][i] - x_t)
        lasso_term = torch.abs(
            [param for name, param in self.G.named_parameters() if name == "weight"][0]
        ).sum()
        self.log("mse", loss, logger=False)
        loss += self.regularization_parameter*lasso_term
        self.log("training_loss", loss, logger=False)
        self.log("lasso_term", lasso_term, logger=False)
        return loss

    def configure_optimizers(self):
        return torch.optim.Adam(self.parameters(), lr=self.learning_rate)

