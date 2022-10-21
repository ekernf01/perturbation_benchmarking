## INSERT_NAME: Modular and customizable regulatory network inference software

This software allows flexible combination of numerous individual components that are common in GRN inference:

- Supervised ML with TF-related features and target-related labels
- Prior and/or learned network structure specifying targets of each TF
- Sharing of information across cell types
- Feature construction using motif analysis or predetermined targets to measure TF activity as opposed to expression

#### Usage

    import whateveriwannanameit

    # Training data must meet certain expectations: e.g. has to have a column "perturbation" specifying which gene is perturbed
    # and a column "perturbation_type" that is "overexpression" or "knockdown" or "knockout".
    report = whateveriwannanameit.validate_training_data(train) 
    
    # Input: adata and (optional) sparse network structure
    grn = whateveriwannanameit.GRN(train, network, ppi) # last 2 args are optional
    
    # Specify a fitting method
    grn.fit(
        method = "linear", 
        confounders = ["donor", "batch"],  # any column in .obs is allowed
        cell_type_sharing = "distinct",    # or "identical" or "similar"
        network_prior = "restrictive",     # or "adaptive" or "informative"
        feature_construction = "identity", # or "activity" or "ppi"
        projection = "factor_graph",       # or "PCA_before_training" or "PCA_after_training"
    )

    # anndata.AnnData is returned
    predictions = grn.predict({("POU5F1", 0), ("NANOG", 0)})