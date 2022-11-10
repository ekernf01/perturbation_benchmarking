## Modular and customizable regulatory network inference software

This software allows flexible combination of numerous individual components that are common in GRN inference:

- Supervised ML with TF-related features and target-related labels
- Network structure specifying targets of each TF (user-provided or learned *de novo*)
- Sharing of information across cell types
- Feature construction using motif analysis or predetermined targets to measure TF activity

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

A simple starting point is to run ridge regression with TF expression as predictors and target expression as the labels. Many other supervised learning methods could form a drop-in replacement for linear models; for example, related literature uses neural networks (DCD-FG), kernel regression (Dynamo), or boosted trees (GRNBOOST). 

- $Y$: target expression
- $X$: TF expression
- $Y \approx X\beta$

It's unclear how to choose the strength of regularization. We want parameters that give optimal estimates of individual causal effects. It seems like MSE-optimal shrinkage would be a sensible choice (e.g. Theorem 1.2 in [Wessel N. van Wieringen's notes](https://arxiv.org/pdf/1509.09169.pdf)), but it is unclear how to estimate this. For now, I'll settle for fast LOOCV (aka PRESS), available via [sklearn's RidgeCV](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.RidgeCV.html). 

#### Sparsity (user-defined)

As in CellOracle, a sparse model can be fit where instead of using all TF's to predict each target, prior knowledge dictates a small number of known TF's to include in the model for each target. This combines easily with supervised learning. If $reg(Y)$ contains user-provided regulators of target gene $Y$, we can use  supervised learning methods (for now, ridge regression) to learn $f$ in $Y \approx f(X_{reg(Y)})$.

#### Sparsity (data-driven)

With or without user-defined network structure, additional pruning can be done to zero out small effects, e.g. with LASSO. The strategy in CellOracle is to keep only the largest RR coefficients (user-defined, but usually about 1e4 total). The model is fit once, then coeffs are selected, and then the model is fit once again on fewer features. 

There are a lot of ways that data-driven and user-defined feature selection methods could be combined. But as of 2022 Nov 10, the only option implemented is:

- Fit ridge regression $Y \approx X_{reg(Y)}\beta$ where $reg(Y)$ contains either all TF's or user-provided regulators of target gene $Y$
- Prune, retaining the largest $N$ coeffs across the whole network, where $N$ is user-provided.

#### Cell type specificity

Some methods use the same model for all training data, but others (e.g. CellOracle) use separate models for each cell type, and some methods specifically focus on adapting to the the data in terms of how similar network structure is across cell types (e.g. csnets, GNAT). As of 2022 Nov 10, we implement two options:

- one model for all cell types
- separate models for each cell type
- (We intend to implement partial sharing of information later.)

There is some unavoidable complexity in how this interacts with sparse network structure, because user-provided or learned structure could be cell-type-specific. For example, CellOracle does separate pruning within each cell type, yielding cell-type-specific network structures even from generic user input. If user-provided structures are cell-type-specific, the cell types would need to be matched between the input network and the input expression data.

Currently, our software *does not* allow user-provided network structure to be cell-type-specific. It *does* allow for cell-type-specific model-fitting and pruning.

#### Feature extraction

In place of TF mRNA levels, some works (notably ARMADA) use genome-wide motif activity.

#### Low-rank structure

Some methods project into a principal subspace to interpret output (e.g. CellOracle). Some methods project into a principal subspace for essentially all modeling (PRESCIENT). Some methods have more bespoke low-rank structure (DCD-FG, sctenifoldknk). In this project it would be useful to offer the following options:

- project into PCA subspace prior to training. "Network" connects PC's to other PC's, not TF's to targets. 
- project targets, but not TFs/features, into PCA subspace prior to training. So network connects TF's to PC's.  
- project into PCA subspace after prediction.

#### Dynamics
#### Interventions

