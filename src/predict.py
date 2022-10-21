import anndata
class GRN:
    """
    Flexible inference of gene regulatory network models.
    """
    def __init__(self, train: anndata.AnnData, networks: dict = None):
        self.train = train 
        self.networks = networks 


    def fit(
        method: str, 
        confounders: list, 
        cell_type_sharing: str = "distinct",   
        network_prior: str = "restrictive",    
        feature_construction: str = "identity",
        projection: str = "factor_graph",      
    ):
        """Fit the model.

        Args:
            method (str): _description_
            confounders (list): _description_
            cell_type_sharing (str, optional): _description_. Defaults to "distinct".
            network_prior (str, optional): _description_. Defaults to "restrictive".
            feature_construction (str, optional): _description_. Defaults to "identity".
            projection (str, optional): _description_. Defaults to "factor_graph".
        """