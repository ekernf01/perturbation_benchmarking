## How to use, re-use, and re-purpose our benchmarking infrastructure

We hope this documentation makes it easy for you to achieve any of the specific goals outlined below. If there's something missing or incorrect, please feel free to file a github issue; we're excited to support anyone repurposing this work. 

### Prereqs and general procedure

These how-to's assume you have already followed our installation instructions. If you have not, please do that first.

Once the software is installed, you will almost always need to: 

- Make a new folder inside `Experiments` named after your experiment
- Populate it with a json file `metadata.json`. We give examples below for specific use-cases.
- Run it in the ggrn conda environment using `do_one_experiment.py`.
- Find results in `Experiments/my_experiment/outputs/evaluationPerPert.parquet`.

##### Code for the last two steps

```bash
cd perturbation_benchmarking
conda activate ggrn
python do_one_experiment.py --experiment_name my_experiment --amount_to_do missing_models --save_trainset_predictions \
    > experiments/my_experiment/stdout.txt 2> experiments/my_experiment/err.txt
```

```python
import pandas as pd
experimental_conditions = pd.read_csv('Experiments/my_experiment/outputs/conditions.csv')
benchmark_results = pd.read_parquet('Experiments/my_experiment/outputs/evaluationPerPert.parquet', engine='pyarrow')
```

##### Expected result 

- `experimental_conditions` should be a small dataframe with one row per experimental condition.
- `benchmark_results` should be a large dataframe with one row per pair (perturbation, experimental condition). The columns will contain performance metrics like mae and metadata. Consult `docs/reference.md` for comprehensive information.

### How to repeat our regression method experiments

Our experiments on different regression methods all start with `1.0_`. We would typically run them like this.

```bash
cd perturbation_benchmarking
conda activate ggrn
for experiment in `ls -1 experiments | grep '1.0_'"`

do
    echo "Starting ${experiment}"
    python do_one_experiment.py --experiment_name $experiment --amount_to_do missing_models \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt
done
```


### How to evaluate a new method

- Make a docker image to run your new method. We have a [separate guide for this](https://github.com/ekernf01/ggrn/tree/main/ggrn_docker_backend).
- Follow the general steps above using the [metadata for our docker demo experiment](https://github.com/ekernf01/perturbation_benchmarking/blob/main/experiments/ggrn_docker_backend/metadata.json) as a starting point.

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

Follow the general procedure discussed above using Experiment `1.1.1_1` (metadata copied below) as an example. The crucial items are `kwargs` and `kwargs_to_expand`.

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

Follow the general procedure discussed above using Experiment `1.8.4_0` (metadata copied below) as a starting point. The crucial items here are `type_of_split` and `data_split_seed`.

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

If you use the `custom` data split, it should be in `Experiments/my_new_experiment/custom_test_sets/<data_split_seed>.json` and it should be formatted as a `json` list containing names of observations to reserve for the test set. In our unit tests, we generate a custom test set file from an AnnData object like this.

```python
import json
custom_test_set = set(adata.obs_names[0:1000])
os.makedirs("custom_test_sets", exist_ok=True)
with open("custom_test_sets/0.json", "w") as f:
    json.dump(list(custom_test_set), f)
```

### How to add a new dataset

See the perturbation data [repo](https://github.com/ekernf01/perturbation_data).

### How to add a new metric

To add your own evaluation metrics, you will need to make a fork of the `perturbation_benchmarking_package` [repo](https://github.com/ekernf01/perturbation_benchmarking_package), edit `evaluator.py`, and install your version prior to running your experiments. 

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

You can add any function by following the same format you see. Results will be included in a column named after the key you add to the dictionary. For example, you could modify the "spearman" item to assess the spearman correlation of predicted vs observed expression (instead of the fold change).
 
```python
{
    "expression_correlation": lambda predicted, observed, baseline: [x for x in spearmanr(observed, predicted)][0],
}
```

##### Installing your version

Navigate to the same folder you ran `git clone` from. Run `conda activate ggrn` and `pip install -e perturbation_benchmarking_experiments`. Then run python and look for your new metric in METRICS. You can test any metric using code similar to this.

```python
from perturbation_benchmarking_package import evaluator
import numpy as np
evaluator.METRICS["mae"](np.array([1,2,3]), np.array([4,5,6]), np.array([7,8,9]))
```

### How to repeat all of our experiments

Our experiments can be run via `source run_experiments.sh &`. Progress can be monitored by inspecting `stdout.txt` and `err.txt` in each experiment's folder. Once the experiments are done, figures can be produced using the R scripts in `make_figures`. 

You are likely to encounter some difficulties.

- Experiments could take a long time (weeks on a typical laptop). We ran experiments bit by bit over a long period, and they are not currently set up to be dispatched in a massively parallel way. 
- The repo is under active development as of December 2023 and may not be entirely stable or may not exactly reproduce our preprint. A list of commit hashes used for version one of our preprint can be found in the `environment` folder, and we plan to make code releases for future preprint versions or journal submissions.
- Making figures requires some common basic R packages like ggplot2 that are not included in our environment setup. Let us know if you have trouble installing them.
