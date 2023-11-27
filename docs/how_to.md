## How to use, re-use, and re-purpose our benchmarking infrastructure

We hope this documentation makes it easy for you to achieve any of the specific goals outlined below. If there's something missing or incorrect, please feel free to file a github issue; we're excited to support anyone repurposing this work. 

### Prereqs

These how-to's assume you have already followed our installation instructions. 

### General procedure

In each case, you will need to 

- Make a new folder named after your experiment.
- Populate it with a `metadata.json`. We give examples below for specific use-cases.
- Run it in the ggrn conda environment.
- Find results in `Experiments/my_experiment/outputs/evaluationPerPert.parquet`.

Code: 

```bash
cd perturbation_benchmarking
conda activate ggrn
python do_one_experiment.py --experiment_name my_experiment --amount_to_do missing_models --save_trainset_predictions \
    > experiments/my_experiment/stdout.txt 2> experiments/my_experiment/err.txt
```

```python
import pandas as pd
benchmark_results = pd.read_parquet('Experiments/my_experiment/outputs/evaluationPerPert.parquet', engine='pyarrow')
```

### How to repeat our regression method experiments

Our experiments on different regression methods all start with `1.0_`. We would typically run them like this.

```bash
cd perturbation_benchmarking
conda activate ggrn
for experiment in `ls -1 experiments | grep '1.0_'"`

do
    echo "Starting ${experiment}"
    echo "Monitor progress:
less experiments/${experiment}/err.txt
less experiments/${experiment}/stdout.txt
"
    python do_one_experiment.py --experiment_name $experiment --amount_to_do missing_models \
        > experiments/$experiment/stdout.txt 2> experiments/$experiment/err.txt
done
```

### How to evaluate a new method

- Make a docker image to run your new method. We have a [separate guide for this](https://github.com/ekernf01/ggrn/tree/main/ggrn_docker_backend).
- Copy and modify the [metadata for our docker demo experiment](https://github.com/ekernf01/perturbation_benchmarking/blob/main/experiments/ggrn_docker_backend/metadata.json).

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

Use Experiment `1.1.1_1` as an example. The crucial items are `kwargs` and `kwargs_to_expand`.

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

Use Experiment `1.8.4_0` as an example. The crucial items here are `type_of_split` and `data_split_seed`.

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