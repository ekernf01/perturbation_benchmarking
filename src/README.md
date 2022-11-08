## Modular and customizable regulatory network inference software

This software allows flexible combination of numerous individual components that are common in GRN inference:

- Arbitrary supervised ML with TF-related features and target-related labels
- Prior and/or learned network structure specifying targets of each TF
- Sharing of information across cell types
- Feature construction using motif analysis or predetermined targets to measure TF activity as opposed to expression

It is a work in progress.

### How to use it

#### Installation

WRITEME

#### Usage

    import predict

    # Training data must meet certain expectations: e.g. has to have a column "perturbation" specifying which gene is perturbed
    # and a column "perturbation_type" that is "overexpression" or "knockdown" or "knockout".
    report = predict.validate_training_data(train) 
    
    # Input: adata and (optional) sparse network structure
    grn = predict.GRN(train, network, ppi) # last 2 args are optional
    
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

### Methods description

This section will describe various options for GRN modeling, along with some possible combinations.

#### Linear modeling

A simple starting point is to run ridge regression with TF expression as predictors and target expression as the labels. Many other supervised learning methods could form a drop-in replacement for linear models; for example, related literature uses GAMs (DCD-FG) or kernel regression (Dynamo). 

- $Y$: target expression
- $X$: TF expression
- $Y \approx X\beta$

It's unclear how to choose the strength of regularization. We want parameters that give optimal estimates of individual causal effects. It seems like MSE-optimal shrinkage would be a sensible choice (e.g. Theorem 1.2 in [Wessel N. van Wieringen's notes](https://arxiv.org/pdf/1509.09169.pdf)), but it is unclear how to estimate this. For now, I'll settle for fast LOOCV (aka PRESS), available via [sklearn's RidgeCV](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.RidgeCV.html). 

#### Sparsity (user-defined)

As in CellOracle, a sparse model can be fit where instead of using all TF's to predict each target, prior knowledge dictates a small number of known TF's to include in the model for each target. This combines easily with supervised learning.

- $reg(Y)$: known regulators of target gene $Y$
- $Y \approx X_{reg(Y)}\beta$

#### Sparsity (data-driven)

With or without user-defined network structure, additional pruning can be done to zero out small coefficients, e.g. with LASSO. The strategy in CellOracle is to keep only the largest coefficients (user-defined, but usually about 1e4 total). The model is fit once, then coeffs are selected, and then the model is fit once again on fewer features. This is compatible with any supervised learning method that produces a feature important measure as part of its output. Kernel regression could be difficult and random forests would require a GENIE3-like importance measure, but ridge regression works fine. Something reasonable could be extracted from DCD-FG even if it's not clear exactly what. 

#### Feature extraction
In place of TF mRNA levels, some works (notably ARMADA) use genome-wide motif activity.

#### Low-rank structure

Some methods project into a principal subspace to interpret output (e.g. CellOracle). Some methods project into a principal subspace for essentially all modeling (PRESCIENT). Some methods have more bespoke low-rank structure (DCD-FG, sctenifoldknk). In this project it would be useful to offer the following options:

- project into PCA subspace prior to training. "Network" connects PC's to other PC's, not TF's to targets. 
- project targets, but not TFs/features, into PCA subspace prior to training. So network connects TF's to PC's.  
- project into PCA subspace after prediction.

#### Dynamics
#### Interventions
#### Cell type specificity

