## GGRN: benchmarking across a grammar of dynamic gene regulatory network inference methods

We describe a "grammar" capable of closely mimicking the predictive methods underlying:

- several NOTEARS variants, including DCDFG
- PRESCIENT
- CellOracle and Dictys
- GENIE3, DynGENIE, and SCENIC
- Dynamo
- ScanBMA
- ARMADA
- LLC (Hyttinen et al 2012)

Our objective is to enable easy combination of different methodological features for benchmarking purposes.

#### Notation 

Let $X_i(P, T)$ be RNA abundance in sample $i$, at time $T$, where genes in the set $P$ were overexpressed continuously starting at $T=0$. A natural approach is to fit $F$ such that $X_i(P, T) - F(X_i(P, T-\delta T))$. 

#### Matching

A prominent obstacle is that RNA-seq destroys samples: $X_i(P, T)$ and $X_i(P, T-\delta T)$ cannot both be observed. The most common work-around is to assume that samples have reached a steady state and to fit $X_i(P, T) \approx F(X_i(P, T))$; this is done by LLC, NOTEARS variants, CellOracle, Dictys, and GENIE3. If time-series data are available, an alternative is to match each sample with a counterpart at the previous time-step, which is simple enough with small numbers of bulk RNA samples (ScanBMA, ARMADA, DynGENIE). With single-cell RNA data, this type of matching is harder; it requires lineage tracing predictions. In one method, PRESCIENT, these are solved by optimal transport matching. Another way to obtain paired measurements across two timepoints is RNA velocity, which is used by Dynamo.

#### Predictive timescale and acyclic penalties

To predict the effect of perturbing gene $j$ to level $z$ (e.g. $z=0$ for a knockout), an obvious choice is to start with control expression $X$, set $X_j \leftarrow z$, and then set $X \leftarrow F(X)$. For models like Dynamo, PRESCIENT, ARMADA, and ScanBMA that are coupled to an explicit time-scale, an entire trajectory can be computed by iterating this process. For steady-state models, the amount of time simulated by a single iteration is unclear. CellOracle leaves this to the user, recommending 3 iterations. In DCD-FG, the handling of this decision is [under discussion](https://github.com/genentech/dcdfg/issues) as of 2022 December 16, but the training assumptions seem to imply iteration until a fixed point is reached. To guarantee existence of a fixed point, training of DCD-FG and related methods includes a differentiable penalty term that is constrained to be 0, purging cycles from the network structure. LLC does not include this type of penalty, but assumes a related "weak stability" condition that rules out feedback loop explosion.  

#### Perturbation persistence 

In PRESCIENT, $X_j=z$ is done once, and only $X = F(X)$ is iterated, corresponding to a transient perturbation (e.g. inducible knockdown). In CellOracle, each iteration includes both $X_j=z$ and $X = F(X)$, corresponding to a persistent perturbation (e.g. knockout). In Dynamo, perturbations also persist; they are incorporated into the entire vector field of derivatives, and a path is selected to minimize the "action", a measure of discordance between the curvature of the path and the predicted vector field of derivatives. LLC also assumes persistent perturbations.

#### TF activity 

Many transcription factors are post-transcriptionally regulated, and the mRNA level may not be the best proxy for TF activity. One method, ARMADA, augments each log-normalized expression profile $X$ with additional coordinates $A$ representing TF activity. $A$ is estimated by fitting $X_{ij} \approx \sum_m N_{jm}A_{mi}$, where $N_{jm}$ counts how many times motif $m$ occurs near the promoter of gene $j$. $N_{jm}$ is constructed from promoter annotations, motif databases, and DNA sequence; it does not depend on expression data except through possibly-novel promoter annotations. Each promoter's expression and motif counts are centered on 0 prior to fitting. 

The idea to quantify TF activity through downstream effects, rather than relying on TF mRNA expression, elegantly accounts for post-transcriptional regulation of TF activity. Mathematically, this can be phrased in terms of low-rank structure, which is described below in the context of other methods. 

#### Regression

The matching methods described above yield a set of paired measurements $\{X_i, Y_i\}$. Every method described by this "grammar" uses supervised ML to fit $F$ such that $Y_i \approx F(X_i)$. Different families are assumed: CellOracle, ScanBMA, Dictys, ARMADA, LLC, NOTEARS, and NOTEARS-LR use linear regression; DCD-FG uses multilayer perceptrons; NOBEARS uses polynomial regression; GENIE uses random forests; Dynamo uses kernel ridge regression; and PRESCIENT uses a neural network. Sometimes these are not trained directly on expression data; they may also involve projection into lower-dimensional latent spaces. 

#### Low dimensional structure

Some methods' estimates of $F$ are heavily affected by assumptions about low-rank structure. Let $Q(X)$ project $X$ into a $D$-dimensional space (via PCA unless otherwise noted). CellOracle requires projection of all output into low-dimensional spaces for biological interpretation of predictions, effectively assuming $F(X) = Q^{-1}(Q(G(X)))$ where $G$ is the Bayesian ridge regression model mentioned above. PRESCIENT and Dynamo estimate derivatives after projection into a 30- or 50-dimensional principal subspace, assuming $F(X) = Q^{-1}(G(Q(X)))$. DCD-FG and NOTEARS-LR each jointly learn an "encoder" $Q$ and a "decoder" $G$ such that $F(X) = G(Q(X))$. For NOTEARS-LR, $G$ and $Q$ are linear, and for DCD-FG, $G$ is linear but $Q$ consists of $D$ separate multilayer perceptrons. NOTEARS-LR and DCD-FG do not use PCA; they learn the encoder and the decoder in a supervised way by backpropagating errors. ARMADA assumes $F(X) = Q^{-1}(G(Q(X)))$, but the "encoder" $Q$ is not learned via PCA or backprop. Rather, it is fixed *a priori* based on motif analysis; in the notation above, $Q(X) = N^{\dagger}X$ where $N^{\dagger}$ is the Moore-Penrose pseudoinverse of the motif-by-promoter count matrix $N$.
 
#### Sparse structure

Estimates of $F$ are also heavily affected by assumptions about sparsity. DCD-FG and NOTEARS variants use an L1 penalty so that each coordinate of $F(X)$ (or $Q(X)$) is a function of just a few coordinates of $X$. ScanBMA uses Bayesian model scoring tools toward the same end, selecting a model with high BIC and averaging its predictions with nearby models above a certain performance threshold. CellOracle and SCENIC use ML methods that yield dense solutions, but they each include a step to prune weak coefficients and re-fit each regression with fewer inputs. LLC uses an L1 penalty to induce sparsity for some applications. PRESCIENT and Dynamo do not assume any sparsity. ARMADA does not assume any additional sparsity beyond what is implied by motif analysis (discussed below).

#### Known or suspected interactions

The specific pattern of connections in the network can be informed by prior knowledge about regulatory relationships -- most often, motif analysis. CellOracle, ARMADA, and SCENIC+ allow regulation of gene $j$ by regulator $k$ only if a matching motif is found near the promoter of $j$ or in a paired enhancer. Earlier iterations of ScanBMA used a similar hard *a priori* threshold, while later ScanBMA papers include the same information as an informative prior on network structure. 

CellOracle, Dictys, and SCENIC+ each include analysis of ATAC-seq data to find motifs and pair enhancers with genes. This may fall outside the scope of our project, but is consistently helpful when comparing learned network structure to ChIP-seq data.

DCD-FG and NOTEARS variants do not use motif analysis. 

#### Cell-type specificity

Many of these methods are demonstrated on fairly homogenous systems like yeast (ScanBMA), hematopoeisis (PRESCIENT, Dynamo) or THP-1 or HUVEC cells (ARMADA). Throughout analysis of each dataset, they use the same models. CellOracle, though, opts to train a separate $F$ on each of many cell types or developmental stages. For CellOracle, the regression, pruning, and refitting are done separately within each cell type. An intermediate option would be to fit cell type-specific models, but shrink each model towards shared structures or effect sizes.

#### Interpretation of randomness

RNA-seq, and especially single-cell RNA-seq, produce data that are noisy. Transcription itself is also stochastic. Most of the methods described here include random elements, which could be interpreted as biological reality or measurement error. The intepretation is not always obvious, but we attempt to faithfully summarize the intent behind each method. CellOracle performs smoothing prior to model fitting and includes no noise in simulations, meaning randomness is interpreted as measurement error. PRESCIENT, Dictys, and DCD-FG each add noise continously in forward simulations, meaning randomness is interpreted as biological reality.

#### Summary

Each method can be described as a combination of: 

- A matching scheme
- A prediction timescale (which may necessitate an acyclic penalty)
- The perturbation persistence
- The regression method
- The low-dimensional structure
- A method of encouraging sparsity
- A method of integrating user-provided known or suspected interactions into the sparsity pattern
- A choice to fit cell type-specific or shared models
- A choice to treat random model components as measurement error or biological reality

Below, we attempt to compress each method into json format, reserving the value `none` for when the original method does not include a given feature.

CellOracle:

    {
        "matching": "steady_state",
        "acyclic_penalty": none,
        "prediction_timescale": "3",
        "perturbation_is_persistent": true,
        "regression_method": "BayesianRidge",
        "low_dimensional_structure": "QiQG",
        "low_dimensional_training": "PCA",
        "sparsity": "prune_refit",
        "prior": "hard_threshold",
        "do_cell_type_specific": true,
        "is_noise_biological": false, 
    }

DCD-FG:

    {
        "matching": "steady_state",
        "acyclic_penalty": "spectral_radius",
        "prediction_timescale": "Inf",
        "perturbation_is_persistent": true,
        "regression_method": "MultiLayerPerceptron",
        "low_dimensional_structure": "GQ",
        "low_dimensional_training": "supervised",
        "sparsity": "L1",
        "prior": "none",
        "do_cell_type_specific": false,
        "is_noise_biological": true, 
    }

PRESCIENT:

    {
        "matching": "optimal_transport",
        "acyclic_penalty": none,
        "prediction_timescale": "real",
        "perturbation_is_persistent": false,
        "regression_method": "NeuralNetwork",
        "low_dimensional_structure": "QiGQ",
        "low_dimensional_training": "PCA",
        "sparsity": "none",
        "prior": "none",
        "do_cell_type_specific": false,
        "is_noise_biological": true, 
    }

GENIE3:

    {
        "matching": "steady_state",
        "acyclic_penalty": none,
        "prediction_timescale": none,
        "perturbation_is_persistent": none,
        "regression_method": "RandomForest",
        "low_dimensional_structure": none,
        "low_dimensional_training": none,
        "sparsity": "none",
        "prior": "none",
        "do_cell_type_specific": false,
        "is_noise_biological": none, 
    }
    
ARMADA:

    {
        "matching": "timeseries",
        "acyclic_penalty": none,
        "prediction_timescale": "real",
        "perturbation_is_persistent": none,
        "regression_method": "linear",
        "low_dimensional_structure": "QiGQ",
        "low_dimensional_training": "fixed",
        "sparsity": "none",
        "prior": "hard_threshold",
        "do_cell_type_specific": false,
        "is_noise_biological": none, 
    }

ScanBMA:

    {
        "matching": "timeseries",
        "acyclic_penalty": none,
        "prediction_timescale": "real",
        "perturbation_is_persistent": none,
        "regression_method": "linear",
        "low_dimensional_structure": none,
        "low_dimensional_training": none,
        "sparsity": "BMA",
        "prior": "probabilistic",
        "do_cell_type_specific": false,
        "is_noise_biological": true, 
    }

#### Limitations

In no particular order:

- This grammar is discrete-time, and can only approximate continuous-time methods (PRESCIENT, Dictys, Dynamo). 
- This grammar cannot specify the exact architectures used in the DCD-FG's factor graph or PRESCIENT's potential function. 
- This grammar cannot capture how PRESCIENT's optimal transport feature is fully, beautifully integrated: PRESCIENT selects $j$ given $i$ not by matching $X_i(T-\delta T)$ to $X_j(T)$, but rather by matching $F(X_i(T-\delta T))$ to $X_j(T)$. This requires jointly optimizing over both $F$ and the matching method's output.
- This grammar also does not describe any data preprocessing, rather assuming log normalized gene expression as input.
- This grammar does not describe any motif analysis or chromatin data analysis. CellOracle, SCENIC+, and Dictys include steps for data-driven pairing of enhancers with promoters. ARMADA emphasizes promoter annotation, models of motif evolution, and motif positioning relative to the promoter. 

#### Implementation

See `GGRN_specs.md`.