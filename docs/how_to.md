## How to use, re-use, and re-purpose our benchmarking infrastructure

We hope this documentation makes it easy for you to achieve any of the goals outlined below. If there's something missing or incorrect, please file a github issue; we're excited to support anyone repurposing this work. 

### Prereqs and general procedure

**Warning!** These how-to's assume you have already followed our [installation instructions](https://github.com/ekernf01/perturbation_benchmarking/blob/main/environment/install.md). Otherwise, you will not be able to follow along and test the commands below.

Once the software is installed, you will almost always need to: 

- Make a new folder inside `experiments` named after your experiment
- Populate it with a json file `metadata.json`. We give examples below for specific use-cases.
- Run it in the ggrn conda environment using the `pereggrn` CLI.
- Find results in `Experiments/my_experiment/outputs/evaluationPerPert.parquet`.

##### Code for the last two steps

```bash
cd perturbation_benchmarking
conda activate ggrn
pereggrn --experiment_name my_experiment --amount_to_do missing_models --save_trainset_predictions \
    > experiments/my_experiment/stdout.txt 2> experiments/my_experiment/err.txt
```

```python
import pandas as pd
import scanpy as sc
experimental_conditions = pd.read_csv('experiments/my_experiment/outputs/conditions.csv')
benchmark_results = pd.read_parquet('experiments/my_experiment/outputs/evaluationPerPert.parquet', engine='pyarrow')
predictions = sc.read_h5ad('experiments/my_experiment/outputs/predictions/0.h5ad')
```

##### Expected result 

- `pereggrn` should populate `experiments/my_experiment/outputs` with gene expression predictions and evaluation results. 
- `experimental_conditions` should be a small dataframe with one row per experimental condition.
- `benchmark_results` should be a large dataframe with one row per pair (perturbation, experimental condition). The columns will contain performance metrics like mae and metadata. Consult `docs/reference.md` for comprehensive information.
- `predictions` should be an AnnData object containing predicted expression. This will change in size depending on the data split. For most splits this is the same size as the test data. For "timeseries", though, it will contain separate predictions for each combination of perturbation, cell type, starting timepoint, and prediction timescale (number of iterations).

##### Troubleshooting

The default behavior of `pereggrn` is optimized for performance, but several non-default options can help a lot with initial use and debugging. Run `pereggrn -h` for usage instructions. In brief:

- If `conda` is not available, try the usual `source "${HOME}/mambaforge/etc/profile.d/conda.sh"` (unix only).
- You can get error tracebacks by not skipping individual bad runs and not using Joblib parallelization. 
- You can save models for later inspection using `save_models`.
- If you are afraid of overfitting, you can monitor train and test-set performance separately using `save_trainset_predictions`.

### How to repeat our network comparison experiments

Our experiments on different regression methods all start with `1.4.3_`. We would typically run them like this.

```bash
cd perturbation_benchmarking
conda activate ggrn
for experiment in `ls -1 experiments | grep '1.4.3_'"`

do
    echo "Starting ${experiment}"
    pereggrn --experiment_name $experiment --amount_to_do missing_models \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt
done
```

The metadata looks like this, and the networks are specified towards the bottom. 

```
{
    "unique_id": "1.4.3_1",
    "nickname": "base_network",
    "readme": "people have published big lists of TF-target or gene-gene relationships, often for GWAS interpretation or reprogramming. Existing benchmarks have limited information content and seldom compare these published network structures directly without introducing confounding factors. For instance, one might ask whether the networks used by CellNet, Mogrify, Irene, and CellOracle are of comparable value in predicting perturbation outcomes. Those methods have been compared, but they each involve many other components that may also affect the outcome, confounding the effect of network structure. This experiment benchmarks many networks using otherwise-equivalent methods to see how much each network helps predict held-out perturbations.",
    "question": "1.4.3",
    "is_active": true,
    "factor_varied": "network_datasets",
    "data_split_seed": [0],
    "color_by": "type_of_split",
    "type_of_split": ["interventional"],
    "facet_by": null,
    "merge_replicates": true,
    "regression_method": "RidgeCV",
    "perturbation_dataset": "nakatake",
    "eligible_regulators": "human_tfs",
    "num_genes": 10000,
    "network_prior": "restrictive",
    "network_datasets": {
        "empty":     { "do_aggregate_subnets": true }, 
        "dense":     { "do_aggregate_subnets": true },    
        "celloracle_human":      { "do_aggregate_subnets": true },
        "gtex_rna":              { "do_aggregate_subnets": true },
        "magnum_compendium_32":  { "do_aggregate_subnets": true },   
        "magnum_compendium_ppi": { "do_aggregate_subnets": true },
        "cellnet_human_Hg1332":  { "do_aggregate_subnets": true },
        "cellnet_human_Hugene":  { "do_aggregate_subnets": true },
        "MARA_FANTOM4":          { "do_aggregate_subnets": true },
        "STRING":                { "do_aggregate_subnets": true },
        "ANANSE_0.5":            { "do_aggregate_subnets": true },
        "ANANSE_tissue_0.5":     { "do_aggregate_subnets": true },
        "humanbase":             { "do_aggregate_subnets": true }
    }
}
```


### How to add a new dataset or network

See the perturbation data [repo](https://github.com/ekernf01/perturbation_data) or network collection [repo](https://github.com/ekernf01/network_collection).

### How to evaluate a network structure with no fancy modeling

Using our infrastructure, you can run an evaluation where instead of training machine learning models, you use networks to predicts positive regulation versus no regulation, similar to the evaluations in [BETS](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008223) figure 6. Follow the general procedure at the top of this file using the metadata below. 

The most important argument is that `regression_method` is set to `"regulon"`. For the `"regulon"` method, no models will be trained. The targets in the provided network will just exactly mirror the fold change of their regulators. For example, in `outputs/predictions/0.h5ad`, we use the `celloracle_human` network, so whenever log-scale FOXA1 expression goes up by 0.5 due to overexpression, expect to see FOXA1's targets in the `celloracle_human` network go up by 0.5 in the predictions. Evaluations compare the predicted targets against the remaining genes, looking for enrichment of changes in the expected direction. In the parquet file with per-perturbation output, look for the evaluation metric called `pvalue_targets_vs_non_targets`, which is the p-value from an ANOVA comparing observed fold change for genes predicted to change (targets) against observed fold change for genes predicted to stay the same (non-targets).

```json
{
    "unique_id": "1.4.4_1",
    "nickname": "network_only",
    "readme": "Are network-connected genes enriched for perturbation responses? This experiment uses network structure alone for prediction, with no training data and all perturbations reserved for evaluation.",
    "question": "1.4.4",
    "factor_varied": "network_datasets",
    "type_of_split": "interventional",
    "desired_heldout_fraction": [1],
    "color_by": null,
    "facet_by": null,
    "regression_method": "regulon",
    "perturbation_dataset": "nakatake",
    "num_genes": 10000,
    "network_datasets": {
        "celloracle_human":      { "do_aggregate_subnets": true },
        "gtex_rna":              { "do_aggregate_subnets": true },
        "magnum_compendium_32":  { "do_aggregate_subnets": true },   
        "magnum_compendium_ppi": { "do_aggregate_subnets": true },
        "cellnet_human_Hg1332":  { "do_aggregate_subnets": true },
        "cellnet_human_Hugene":  { "do_aggregate_subnets": true },
        "MARA_FANTOM4":          { "do_aggregate_subnets": true },
        "STRING":                { "do_aggregate_subnets": true },
        "ANANSE_0.5":            { "do_aggregate_subnets": true },
        "ANANSE_tissue_0.5":     { "do_aggregate_subnets": true },
        "humanbase":             { "do_aggregate_subnets": true }
    }
}
```

### How to evaluate a new method

- Make a docker image to containerize your new method. We have a [separate guide for this](https://github.com/ekernf01/ggrn/tree/main/ggrn_docker_backend).
- Follow the general steps given at the top of this file using the [metadata for our docker demo experiment](https://github.com/ekernf01/perturbation_benchmarking/blob/main/experiments/ggrn_docker_backend/metadata.json) as a starting point.

```json
{
    "unique_id": "ggrn_docker_backend",
    "nickname": "ggrn_docker_backend",
    "readme": "Test of the ggrn backend that runs a user-specified docker container.",
    "question": "None", 
    "is_active": true,
    "factor_varied": "regression_method",  
    "color_by": null,
    "facet_by": null,
    "perturbation_dataset": "nakatake",
    "kwargs": {
        "my_sweepy_hyperparameter": [0, 1],
        "my_constant_hyperparameter": 0
    },
    "kwargs_to_expand": ["my_sweepy_hyperparameter"],
    "regression_method":"docker____ekernf01/ggrn_docker_backend" 
}
```

### How to run a hyperparameter sweep

Follow the general procedure discussed at the top of this file using Experiment `1.1.1_1` (metadata copied below) as an example. The crucial items are `kwargs` and `kwargs_to_expand`; if it is unclear, you can read about them in `docs/reference.md`. You can combine this with the Docker example above. 

```json
{
    "unique_id": "1.1.1_1",
    "nickname": "hyperparam_sweep",
    "readme": "This is a simple sweep over the regularization parameter for LASSO regression.",
    "question": "1.1.1",
    "is_active": true,
    "factor_varied": "alpha",
    "regression_method": "LASSO",
    "kwargs_to_expand": ["alpha"],
    "kwargs":{
        "alpha": [0, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.01, 0.1, 1, 10, 100, 1000]
    },
    "facet_by": null,
    "color_by": null,
    "perturbation_dataset": "nakatake",
    "eligible_regulators": "human_tfs"
}
```

### How to split the data differently

Follow the general procedure discussed at the top of this file using Experiment `1.8.4_0` (metadata copied below) as a starting point. The crucial items here are `type_of_split` and `data_split_seed`; if it is unclear, you can read about them in `docs/reference.md`.

```json
{
    "unique_id": "1.8.4_0",
    "nickname": "data split",
    "readme": "This is intended to test whether the problem actually gets harder when we split the data 'unevenly', with no common perturbations between train and test data.",
    "question": "1.8.4",
    "is_active": true,
    "facet_by": "type_of_split",
    "color_by": "regression_method",
    "factor_varied": "data_split_seed",
    "regression_method": [ "median", "mean", "RidgeCV"],
    "type_of_split": ["simple", "interventional", "stratified"],
    "data_split_seed": [0,1,2],
    "num_genes": 10000,
    "eligible_regulators": "human_tfs",
    "merge_replicates": false,
    "perturbation_dataset": "nakatake",
    "network_datasets": {
        "dense":{}
    }
}
```

If you use the `custom` data split, custom data split info should be in `Experiments/my_new_experiment/custom_test_sets/<data_split_seed>.json` and it should be formatted as a `json` list containing names of observations to reserve for the test set. In our unit tests, we generate a custom test set file from an AnnData object like this.

```python
import json
custom_test_set = set(adata.obs_names[0:1000])
os.makedirs("custom_test_sets", exist_ok=True)
with open("custom_test_sets/0.json", "w") as f:
    json.dump(list(custom_test_set), f)
```

Here is an example of the way results can differ in response to different data splits. 

![A heatmap showing different results from 'interventional', 'simple', and 'stratified' splits](fig_data_splitting.pdf)

This example can be constructed by running experiment `1.8.4_0` and using the following R code for plotting. 

```r
library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
setwd("perturbation_benchmarking/make_figures/")
source("plotting_functions.R") # Defines collect_experiments and heatmap_all_metrics
X = collect_experiments("1.8.4_0") 
X$x = X$regression_method
heatmap_all_metrics(X, facet1 = "data_split_seed", facet2 = "type_of_split", compare_across_rows = FALSE)
ggsave('fig_data_splitting.pdf', width = 8, height = 8)
```


### How to add a new evaluation metric

To add your own evaluation metrics, you will need to make a fork of the [pereggrn repo](https://github.com/ekernf01/perturbation_benchmarking_package), edit `evaluator.py`, and install your version prior to running your experiments. 

##### Making the fork

Use github's [interface](https://docs.github.com/en/get-started/quickstart/fork-a-repo). Use `git clone` to download your fork.

##### Editing the code

Find `evaluator.py`. Near the top is a dict `METRICS` containing our evaluation functions. It should like this.

``` python
METRICS = {
    "mae":                          lambda predicted, observed, baseline: np.abs(observed - predicted).mean(),
    "mse":                          lambda predicted, observed, baseline: np.linalg.norm(observed - predicted)**2,
    "spearman":                     lambda predicted, observed, baseline: [x for x in spearmanr(observed - baseline, predicted - baseline)][0],
    "proportion_correct_direction": lambda predicted, observed, baseline: np.mean(np.sign(observed - baseline) == np.sign(predicted - baseline)),
    "mse_top_20":                   lambda predicted, observed, baseline: mse_top_n(predicted, observed, baseline, n=20),
    "mse_top_100":                  lambda predicted, observed, baseline: mse_top_n(predicted, observed, baseline, n=100),
    "mse_top_200":                  lambda predicted, observed, baseline: mse_top_n(predicted, observed, baseline, n=200),
}
```

The basic inputs are, for a given predicted profile, the predicted expression, the observed expression, and baseline expression (by default, the average of the control samples from the training data). You can add any function by following the same format you see. Results will be included in a column named after the key you add to the dictionary. For example, you could add an "expression_correlation" item to assess the spearman correlation of predicted vs observed expression (as opposed to our existing option "spearman" which compares predicted vs observed fold change over baseline).
 
```python
{
    "expression_correlation": lambda predicted, observed, baseline: [x for x in spearmanr(observed, predicted)][0],
}
```

##### Installing your version

Navigate to the same folder you ran `git clone` from. Run `conda activate ggrn` and `pip install -e perturbation_benchmarking_experiments`. Then run python and confirm that your new metric appears in METRICS. You can manually test any metric using code similar to this, and you can add automated tests to `tests/test_evaluator.py`.

```python
from perturbation_benchmarking_package import evaluator
import numpy as np
evaluator.METRICS["mae"](np.array([1,2,3]), np.array([4,5,6]), np.array([7,8,9]))
```

### How to repeat all of our experiments

Our experiments can be run via `./run_experiments.sh &`. Progress can be monitored by inspecting `stdout.txt` and `err.txt` in each experiment's folder. Once the experiments are done, figures can be produced using the R scripts in `make_figures`. 

You are likely to encounter some difficulties.

- Experiments could take a long time (weeks to months on a typical laptop). We ran experiments bit by bit over a long period, and they are not currently set up to be dispatched to cloud or cluster resources in a massively parallel way. If it's worth the investment to you, a good option might be to convert `run_experiments.sh` into a SnakeMake or WDL pipeline.
- The repo is under active development as of February 2024 and may not be entirely stable or may not exactly reproduce our preprint. A list of commit hashes used for version one of our preprint can be found in the `environment` folder, and we plan to make code releases for future preprint versions or journal submissions. Certain output may never be repeatable: notably, DCD-FG is not deterministic.
- Making figures requires some commonly used R packages like ggplot2 that are not included/versioned in our environment setup. Please let us know if you have trouble installing them.
