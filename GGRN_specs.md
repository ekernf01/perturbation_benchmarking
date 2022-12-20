### GGRN software specifications

From the grammar set out in `GGRN.md`, many of the combinations that can be specified either make no sense, or would be very complex to implement. The current software is has two separate backends with a third planned. Each has a different subset of features.

### Backend 1: basics 

The initial implementation requires the steady-state matching method. It offers flexible regression, rudimentary priors, rudimentary cell type specificity, and predictions after just one step forward. It has no other scheme for matching controls to treatment samples, no acyclic penalty, no low-dimensional structure, and no biological noise. It can be summarized by the following JSON.

    {
        "matching": "steady_state",
        "acyclic_penalty": none,
        "prediction_timescale": "1",
        "perturbation_is_persistent": none,
        "regression_method": ["bunch of sklearn options"],
        "low_dimensional_structure": none,
        "low_dimensional_training": none,
        "sparsity": ["none", "prune_and_refit", "built_in"],
        "prior":    ["none", "hard"],
        "do_cell_type_specific": false,
        "is_noise_biological": false
    }

Initial benchmarks across the above attributes were entirely negative: the mean of the training data outperforms all other methods. This motivated backend 2.

### Backend 2: DCD-FG 

One of the most promising benchmarks in the literature is on DCD-FG, so we are incorporating DCD-FG into GGRN. This offers a thin wrapper with no features beyond the original DCD-FG implementation. Fortunately, the original DCD-FG implementation spans a variety of related methods. 

    {
        "matching": "steady_state",
        "acyclic_penalty": ["spectral", "exp"],
        "prediction_timescale": "Inf",
        "perturbation_is_persistent": true,
        "regression_method": ["linear", "multilayer_perceptron"],
        "low_dimensional_structure": "EncoderDecoder",
        "low_dimensional_training": "supervised",
        "sparsity": ["built_in"],
        "prior":    ["none"],
        "do_cell_type_specific": false,
        "is_noise_biological": true 
    }

Initial results from modifying the DCD-FG repo show that DCD-FG performs worse than a simple Gaussian baseline that uses no information about interventions (and makes the same prediction for every held-out condition). We will continue to work with DCD-FG, hoping for a believable win over baseline, but we also plan to extend GGRN to explore key elements missing from DCD-FG. 

### Backend 3: autoregressive with matching

We next plan to implement the following features.

    {
        "matching": ["random_control", "closest_control"]
        "acyclic_penalty": none,
        "prediction_timescale": "real",
        "perturbation_is_persistent": true,
        "regression_method": ["linear", "multilayer_perceptron"],
        "low_dimensional_structure": ["none", "QiGQ"],
        "low_dimensional_training": ["supervised", "PCA", "fixed"],
        "sparsity": ["built_in"],
        "prior":    ["eric_is_confused"],
        "do_cell_type_specific": false,
        "is_noise_biological": true 
    }

Formally, we will minimize

$$ L(X) = \sum_{i\in treated} ||(R \circ G \circ Q \circ P_i)^S(X_{M(i)}) - X_i||^2 + \\
\sum_{i\in steady} ||(R \circ G \circ Q\circ P_i)(X_{i}) - X_i||^2 + \\ 
J(G, Q, R) $$

where:

- $treated$ is a set of indices for samples that have undergone treatment.
- $steady$ is a set of samples that is assumed to have reached a steady state (usually, all controls).
- $Q$ is a projection matrix and $R$ is some attempt to invert $Q$. $Q$ and $R$ can be learned by PCA, OR learned by backprop, OR specified by the user as e.g. motif counts + pseudoinverse, OR set to the identity matrix. Maybe eventually we could have some rows fixed by the user and others optimized.
- $G$ predicts a single step forward in time by $T/S$ hours, where $T$ is the duration of treatment.
- $P_i$ enforces interventions on genes perturbed in sample $i$.
- $F^S(X) = F(F(...F(X)))$ ($S$ iterations of $F(X)$).
- $M(i)$ is the index of a control sample matched to treated sample $i$. $M(i)$ can be implemented by choosing a random control, OR by choosing the closest control, OR maybe eventually by choosing $M$ to minimize $L(X)$, OR by optimal transport. 
- $J$ is a regularizer, e.g. an L1 penalty on entries of matrices representing $G$ if $G$ is linear.

This framework can be trained on time-series data, interventional data, or a mixture. This will already be a distinctive advantage whether or not it actually wins benchmarks. For time-series data, $S$ should be adjusted per-sample so it is proportional to the time elapsed. 

It is not clear to me whether and how prior information should be included, but the most likely method is to include motif information the way ARMADA does: by fixing $Q$ and $R$. One can imagine various schemes for sharing information across cell types but I think we may omit that from the project scope for now. 

We will probably implement this in Pytorch, so that we can fit basically any differentiable functions for R, P, G, Q, and J. 