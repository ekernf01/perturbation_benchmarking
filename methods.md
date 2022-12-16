## GGRN: benchmarking across a grammar of dynamic gene regulatory network inference methods

We describe a framework capable of closely mimicking the predictive methods underlying:

- several NOTEARS variants, including DCDFG
- PRESCIENT
- CellOracle and Dictys
- GENIE3, DynGENIE, and SCENIC
- Dynamo
- ScanBMA
- ARMADA

Our objective is to enable easy combination of different methodological features for benchmarking purposes.

#### Notation 

Let $X_i(P, T)$ be RNA abundance in sample $i$, at time $T$, where genes in the set $P$ were overexpressed continuously starting at $T=0$. A natural approach is to fit $F$ such that $X_i(P, T) \approx F(X_i(P, T-\delta T))$. 

#### Matching

A prominent obstacle is that RNA-seq destroys samples: $X_i(P, T)$ and $X_i(P, T-\delta T)$ cannot both be observed. The most common work-around is to assume that samples have reached a steady state and to fit $X_i(P, T) \approx F(X_i(P, T))$; this is done by NOTEARS variants, CellOracle, Dictys, and GENIE3. If time-series data are available, an alternative is to match each sample with a counterpart at the previous time-step, which is simple enough with small numbers of bulk RNA samples (ScanBMA, ARMADA, DynGENIE). With single-cell RNA data, this type of matching is harder; it requires lineage tracing predictions. In one method, PRESCIENT, these are solved by optimal transport matching. Another way to obtain paired measurements across two timepoints is RNA velocity, which is used by Dynamo.

#### Predictive timescale and acyclic penalties

To predict the effect of perturbing gene $j$ to level $z$ (e.g. $z=0$ for a knockout), an obvious choice is to start with control expression $X$, set $X_j \leftarrow z$, and then set $X \leftarrow F(X)$. For models like Dynamo, PRESCIENT, ARMADA, and ScanBMA that are coupled to an explicit time-scale, an entire trajectory can be computed by iterating this process. For steady-state models, the amount of time simulated by a single iteration is unclear. CellOracle leaves this to the user, recommending 3 iterations. In DCD-FG, the handling of this decision is [under discussion](https://github.com/genentech/dcdfg/issues) as of 2022 December 16, but the training assumptions seem to imply iteration until a fixed point is reached. To guarantee existence of a fixed point (other than the fixed point `[np.nan for _ in range(G)]`, rofl), training of DCD-FG and related methods includes a penalty term that purges cycles from the network structure.

#### Perturbation persistence 

In PRESCIENT, $X_j=z$ is done once, and only $X = F(X)$ is iterated, corresponding to a transient perturbation. In CellOracle, iterations include both $X_j=z$ and $X = F(X)$. The original ARMADA does not implement this feature, although the authors mention it as a reasonable extension of their work. 

????? Dynamo ?????

#### TF activity 

Many transcription factors are post-transcriptionally regulated, and the mRNA level may not be the best proxy for TF activity. One method, ARMADA, augments each expression profile $X$ with additional coordinates $W$ representing TF activity. $W$ is computed as $XC$, where $C$ is a matrix and $C_{jm}$ counts how many times motif $m$ occurs near gene $j$. ??? NOT SURE ABOUT DETAILS. ???

#### Regression

The matching assumptions described above yield a set of paired measurements $\{X_i, Y_i\}$. Every method described by this "grammar" uses supervised ML to fit $F$ such that $Y_i \approx F(X_i)$. Different families are assumed: CellOracle, ScanBMA, Dictys, ARMADA, NOTEARS, and NOTEARS-LR use linear regression; DCD-FG uses multilayer perceptrons; NOBEARS uses polynomial regression; GENIE uses random forests; Dynamo uses kernel ridge regression; and PRESCIENT uses a neural network. 

#### Low dimensional structure

Estimates of $F$ are heavily affected by assumptions about low-rank structure. Let $Q(X)$ project $X$ into a $D$-dimensional space (via PCA unless otherwise noted). CellOracle requires projection of all output into low-dimensional spaces for biological interpretation of predictions, effectively assuming $F(X) = Q^{-1}(Q(G(X)))$ where $G$ is the supervised ML method mentioned above. PRESCIENT conducts all modeling after projection into a 50-dimensional principal subspace, assuming $F(X) = Q^{-1}(G(Q(X)))$. DCD-FG and NOTEARS-LR each jointly learn an "encoder" $Q$ and a "decoder" $G$ such that $F(X) = G(Q(X))$. For NOTEARS-LR, $G$ and $Q$ are linear, and for DCD-FG, $G$ is linear but $Q$ consists of $D$ separate multilayer perceptrons. NOTEARS-LR and DCD-FG do not use PCA; they learn both the "encoder" and the "decoder" in a supervised way by backpropagating errors. 
 
#### Sparse structure

Estimates of $F$ are also heavily affected by assumptions about sparsity. DCD-FG and NOTEARS variants use an L1 penalty so that each coordinate of $F(X)$ (or $Q(X)$) is a function of just a few coordinates of $X$. CellOracle and SCENIC use ML methods that yield dense solutions, but they each include a step to prune weak coefficients and re-fit each regression with fewer inputs. ScanBMA and ARMADA do ?????. Strangely, PRESCIENT and Dynamo do not assume any sparsity.

#### Prior knowledge

The specific pattern of sparsity can be informed by prior knowledge about regulatory relationships -- most often, motif analysis. CellOracle and SCENIC+ allow regulation of gene $j$ by regulator $k$ only if a matching motif is found in the promoter of $j$ or in a paired cis-regulatory element. Earlier iterations of ScanBMA used a similar hard *a priori* threshold, while later ScanBMA papers include the same information as an informative prior on network structure. 

CellOracle, Dictys, and SCENIC+ each include analysis of ATAC-seq data to find motifs and pair CRE's with genes. This may fall outside the scope of our project, but is consistently helpful when comparing learned network structure to ChIP-seq data.

DCD-FG and NOTEARS variants do not use prior knowledge about structure. 

#### Cell-type specificity

Many of these methods are designed for fairly homogenous systems like yeast (ScanBMA), hematopoeisis (PRESCIENT, Dynamo) or THP-1 differentiation (ARMADA). Throughout analysis of each dataset, they use the same models. Other methods (CellOracle, ????) opt to train a separate $F$ on each of many cell types or developmental stages. For CellOracle, the regression, pruning, and refitting are done separately within each cell type. An intermediate option would be to fit cell type-specific models, but shrink each model towards shared structures or effect sizes.

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
- A source of prior knowledge and a method of integrating it into the sparsity pattern
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
        "low_dimensional_structure": "QGQ",
        "low_dimensional_training": "PCA",
        "sparsity": "prune_refit",
        "prior": "hard",
        "do_cell_type_speciifc": true,
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
        "do_cell_type_speciifc": false,
        "is_noise_biological": true, 
    }

PRESCIENT:

    {
        "matching": "optimal_transport",
        "acyclic_penalty": none,
        "prediction_timescale": "real",
        "perturbation_is_persistent": false,
        "regression_method": "NeuralNetwork",
        "low_dimensional_structure": "QGQ",
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

#### Limitations

This framework is discrete-time, and can only approximate continuous-time methods (PRESCIENT, Dictys, Dynamo). It mostly cannot capture latent variables. It cannot specify the exact architectures used in the DCD-FG's factor graph or PRESCIENT's potential function. It cannot capture how PRESCIENT's optimal transport feature is fully, beautifully integrated: PRESCIENT selects $j$ given $i$ not by matching $X_i(T-\delta T)$ to $X_j(T)$, but rather by matching $F(X_i(T-\delta T))$ to $X_j(T)$, which requires simultaneous optimization of $F$ and the matching output. 

#### Implementation

Many of the combinations that can be specified in this grammar either make no sense, or would be very complex to implement. The initial basic implementation provides only steady-state matching, flexible regression, rudimentary priors, rudimentary cell type specificity, and predictions after just one step forward.  It has no matching scheme, no acyclic penalty, no low-dimensional structure, and no biological noise. It can be summarized by the following not quite machine-readable JSON.

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

Initial benchmarks across the above attributes were entirely negative: the mean of the training data outperforms all other methods. One of the most promising benchmarks in the literature is on DCD-FG, so ongoing work focuses on incorporating DCD-FG into GGRN. This will allow benchmarks across the following attributes. 

    {
        "matching": "steady_state",
        "acyclic_penalty": ["spectral", "exp"],
        "prediction_timescale": "Inf",
        "perturbation_is_persistent": true,
        "regression_method": ["linear", "multilayer_perceptron"],
        "low_dimensional_structure": "QG",
        "low_dimensional_training": "supervised",
        "sparsity": ["built_in"],
        "prior":    ["none"],
        "do_cell_type_specific": false,
        "is_noise_biological": true 
    }

We hope DCD-FG will clearly beat trivial baselines, providing a path forward to explore matching schemes and the use of prior knowledge.