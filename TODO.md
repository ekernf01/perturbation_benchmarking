## To do

- Include some single-gene plots to illustrate spearman = 0.2 vs spearman = 0.3
- Get more test datasets??
- Update unit tests oh gawd
- **Why are some genes easier to predict than others? What makes them easier?**
- **How are these tools used for GWAS or genetics?**
- **How are these tools used for CFE?**

### Current experiments

- **baseNetwork**: people have published big fat lists of TF-target or gene-gene relationships, often for GWAS interpretation or reprogramming. Existing benchmarks have limited information content and seldom compared these published network structures directly without introducing confounding factors. For instance, one might ask whether the networks used by CellNet, Mogrify, Irene, and CellOracle are of comparable value in predicting reprogramming outcomes. Those methods have been compared, but they each involve many other components that may also affect the outcome, confounding the effect of network structure. 
- **cellTypeSpecific**: Do networks inferred for the cell type of interest work better than global networks or networks from the wrong cell types? We try this with CellNet and ANANSE cell-type-specific networks, ESC versus others.
- **clusterSpecific**: CellOracle fits separate linear models on each cluster. Competitors like Dynamo fit a single model across all input data. Which approach works best? How many clusters should be used? The answer may depend on both dataset size and biological diversity, so it'll hard be to extract general lessons, but I still want to know the answer to this. 
- **pruning**: Some models use sparse network structures and some dense. Preliminary experiments indicate that restrictive sparse structures (pre-specified by the user) are not competitive in terms of predicting held-out perturbations. 
- **datasetSize** Does the sparse versus dense tradeoff come out differently on smaller datasets? Alexis posits that sparse models may perform better on smaller networks. Take a good performer from the sparse models and test it against dense on a variety of sizes of training dataset. 

### Future experiments and questions

- A variety of real datasets: CMAP, Fantom4, Ko EST, and maybe one velocity dataset. But why? What's the most interesting thing we learn from this?
- A variety of simulated datasets, simulated from the Ko data with various network structures. With no misspecification, how is the performance?

#### Network structures to include on top of what we have already:

- Something related to IRENE??
- random controls with equal density to each actual network
- optimal structure, see below

#### Optimal Structure Experiment

What's the best possible outcome if we got the network structure from an oracle?
We want min_G(J(X, Y, G)) where J is CellOracle's error from training on X with base network G and testing on Y. 
How might we compute this? Can't try all networks.
Dan Peng offers a clever approximation: take the network structure from the TEST dataset Y. 
But don't take the regression coefficients from Y.
Does it beat everything else? If so, then by how much?

### Perturbations to include

There's a huge universe to explore here. Some items to consider:

- METABRIC (cf Diptavo)
- CD4 (cf Josh W)
- CMAP
- FANTOM4
- CellOracle perturb-seq
- List of papers using CellOracle

