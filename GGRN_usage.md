### How to use it

#### Installation

WRITEME

#### Usage

    import ggrn

    # Training data must meet certain expectations: e.g. has to have a column "perturbation" specifying which gene is perturbed
    # and a column "perturbation_type" that is "overexpression" or "knockdown" or "knockout".
    report = ggrn.validate_training_data(train) 
    
    # Input: adata and (optional) sparse network structure
    grn = ggrn.GRN(train, network, ppi) # last 2 args are optional
    
    # Simple feature construction -- use the mRNA levels of the TF's as a proxy for activity
    grn.extract_features(method = "tf_rna")

    # Specify a fitting method
    grn.fit(
        method = "linear", 
        confounders = ["donor", "batch"],  # any column in .obs is allowed
        cell_type_sharing = "distinct",    # or "identical" or "similar"
        network_prior = "restrictive",     # or "adaptive" or "informative"
        pruning_strategy = "none",         # or "prune_and_refit" or eventually maybe others
        pruning_parameter = None,          # e.g. lasso penalty or total number of nonzero coefficients
        projection = "factor_graph",       # or "PCA_before_training" or "PCA_after_training"
    )

    # anndata.AnnData is returned
    predictions = grn.predict({("POU5F1", 0), ("NANOG", 0)})

