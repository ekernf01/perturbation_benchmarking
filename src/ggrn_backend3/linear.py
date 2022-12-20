import torch
import pytorch_lightning as pl

def P(x: torch.tensor, perturbations: list) -> torch.tensor:
    """Enact a perturbation.
    Args:
        x (torch.tensor): Expression vector
        p (list): list of tuples (target gene index, level), e.g. (5, 0) where gene 5 is Nanog and you knocked out Nanog.

    Returns:
        torch.tensor: Expression with targeted gene set to prescribed level. 
    """
    for p in perturbations:
        x[p[0]] = p[1]
    return x


class LinearAutoregressive(pl.LightningModule):
    def __init__(self, S, G):
        super().__init__()
        self.G = G
        self.S = S
    
    def forward(self, input_batch):
        observations, matched_controls, perturbations, is_treated, is_steady_state = input_batch
        loss = 0
        for i in range(len(x)):
            if is_steady_state[i]:
                # (electric slide voice) one hop this time. *BEWMP*
                x_t = observations[i]
                x_t = torch.matmul(self.G, P(x_t, perturbations[i]))
                loss += torch.linalg.norm(observations[i] - x_t)
            if is_treated[i]:
                # (electric slide voice) s hops this time. *BEWMP* *BEWMP* *BEWMP* *BEWMP* 
                x_t = matched_controls[i]
                for _ in range(self.S-1):
                    x_t = torch.matmul(self.G, P(x_t, perturbations[i]))
                loss += torch.linalg.norm(observations[i] - x_t)
        lasso_term = torch.abs(self.G).sum()
        loss += lasso_term
        return loss



