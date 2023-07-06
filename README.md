## A systematic comparison of computational methods for in silico genetic perturbation

This repo contains benchmark experiments to evaluate various strategies for predicting in detail the outcome of perturbation experiments. To reproduce our results, you'll need to install our code, download our data to the expected relative path, run our Experiments, and finally re-create the figures.

### Related infrastructure

This project is tightly coupled with our collections of data, our GGRN package for network inference, and a companion package containing benchmarking infrastructure. It will not work without all of those components. Our installation script will attempt to set them all up automatically. You can also learn more about them at the following links.

- Perturbation data, the network collection, and some accessory data (e.g. a list of TF's) are on Zenodo with DOI `10.5281/zenodo.8071809`.
    - Our code expects each of those three folders to be unzipped and placed adjacent to this repo.
    - Use our [perturbation loader](https://github.com/ekernf01/load_perturbations) and [network loader](https://github.com/ekernf01/load_networks) to easily access and validate data from Python.
- [GGRN](https://github.com/ekernf01/ggrn) offers flexible combination of different features for regulatory network inference
- Our [perturbation benchmarking package](https://github.com/ekernf01/perturbation_benchmarking_package) helps conduct the `Experiment`s that are specified in this repo.
- Certain additional experiments are implemented in [our fork of DCD-FG](https://github.com/ekernf01/dcdfg).

### Resource requirements and installation

250GB of disk space and 64GB of RAM is enough resources to run everything except certain tree-based models. Certain models nominally require GPU's, but we have been able to run most experiments using a CPU, sometimes by making minimal changes to Pytorch code. See the GGRN repo for details on specific methods.

The project is written in Python. We use Conda + Mamba to manage most dependencies. We offer either a flexible install or an exact install of our environment. Each option has a CPU version and a GPU version, though the GPU version is not well-tested.
 
#### Flexible install

Some flexibility in exact requirements is useful in case of different operating systems. To install roughly the same versions we use, we start with the dependencies given in `environment/conda_inputs.yaml`. We also use a few packages only available via pip. Here is the code.

```bash
mamba env create --name ggrn --file environment/conda_inputs.yaml
conda activate ggrn
pip install vl-convert-python
# Avoid interfering with gimmemotifs 0.17 or with the PyTorch install by using --no-deps
# The deps should be taken care of by the above.
mamba install 'pandas>=1.2' --no-deps
pip install cell-gears==0.0.4  --no-deps
pip install celloracle==0.12.0 --no-deps
pip install prescient==0.1.0   --no-deps 
pip install geomloss==0.2.3    --no-deps 
pip install git+https://huggingface.co/ctheodoris/Geneformer@f0b6641f --no-deps
pip install git+https://github.com/bowang-lab/scFormer@2df344a --no-deps
pip install 'scib>=1.0.3' --no-depsfor p in load_networks load_perturbations ggrn_backend2 ggrn_backend3 ggrn perturbation_benchmarking_package do
    pip install git+https://github.com/ekernf01/${p}
done
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
pip install .
cd ..
```

#### Exact reproduction

Exactly reproducing our environment requires Ubuntu 20.04. The environment can be reproduced via `environment/set_up_environment.sh`. 

### Experiments

The benchmarks in this project are composed of small, structured folders called Experiments. Code and metadata supporting each Experiment resides in a folder in `experiments`, and an individual Experiment can be run by manually editing `run_experiments.sh` or by providing args to `do_one_experiment.py`. The most important parts of each Experiment are:

- `outputs`: these results are auto-generated by our software based on user input and our data collections.
- `metadata.json`: A complete specification of the Experiment, provided by the user. 

Much more detail on metadata and outputs below.

It can be hard to understand how all our experiments relate to one another, software-wise or science-wise. Run `gather_experiment_metadata.py` to produce a quick summary table describing all experiments, or read `all_experiments.tsv` for the most recent such summary. 

#### Outputs

Here is an abridged, annotated description of Experiment outputs.

```bash
├── perturbations # Same as targets but stratified by perturbation instead
├── targets 
│   ├── predictability_vs_in-degree.svg # Display of MAE of groups of targets stratified by in-degree in our networks.
│   ├── variety_in_predictions # Histogram meant to answer, "Are the predictions roughly constant?"
│   ├── enrichr_on_best # Enrichr pathway analysis of best-predicted targets for each condition in this experiment.
│   ├── best # Scatterplots of predicted vs observed for the best-predicted targets.
│   ├── random # Scatterplots of predicted vs observed for a few randomly chosen targets.
│   └── worst # Scatterplots of predicted vs observed for the worst-predicted targets.
├── experiments.csv # All combinations of values provided in the metadata. Would be better named "conditions.csv". 
├── fitted_values # Predictions on training data ...
│   ├── 0.h5ad # ... from row 0 of experiments.csv
│   ├── 1.h5ad # ... from row 1 of experiments.csv
│   ├── ...
├── genes_modeled.csv # The genes included in this experiment.
├── mae.svg # Mean absolute prediction error for each test set observation
├── new_experiments.csv # This is generated and compared to any existing experiments.csv to prevent confusion upon editing metadata.
├── predictions # Predictions on test data 
│   ├── 0.h5ad 
│   ├── 1.h5ad
│   ├── ...
├── evaluationPerPert.parquet # Table of evaluation metrics listed separately for each observation in the test data, readable by e.g. pandas.read_parquet()
├── evaluationPerTarget.parquet # Table of evaluation metrics listed separately for each feature in the test data, readable by e.g. pandas.read_parquet()
└── trainset_performance # Same as parquet files above, but for train-set
    ├── evaluationPerPert.parquet
    ├── evaluationPerTarget.parquet
```

#### Metadata and designing an Experiment

For easy/immediate use, we recommend finding an Experiment similar to what you want to do and then editing it. Metadata will be validated by `perturbation_benchmarking_package.experimenter.validate_metadata()`, and this will provide useful error messages for simple issues such as required args that are missing. For more explicit documentation, read on. 

Experiment metadata files are JSON dictionaries. Most simple entries can be either a single value, or a list. If a list is provided, the experiment is run once for each item in the list. If multiple keys have lists, all combinations will be used. 

With apologies, many metadata keys have idiosyncratic formatting and meaning. 

- `perturbation_dataset` describes a dataset using the same names as our perturbation dataset collection. Only one dataset is allowed per Experiment. 
- `readme` describes the purpose of the experiment. `nickname` conveys the essence curtly. 
- `unique_id` must match the folder the Experiment is in.
- `question` refers to `guiding_questions.txt` in this repo. 
- `is_active` must be `true` or the experiment won't run. 
- `skip_bad_runs`, if `true`, will allow experiments to continue if one condition encounters an error. Set this to `false` for easier debugging.
- `refers_to` points to another Experiment. If A refers to B, then all key/value pairs are copied from B's metadata unless explicitly provided in A's metadata. You may not refer to an experiment that already refers to something. You may not refer to multiple experiments.
- `kwargs` is a dict of keyword args passed on to GEARS, or DCD-FG, or any method [wrapped via Docker](https://github.com/ekernf01/ggrn_docker_backend).
- `baseline_condition` is a number, most often 0. This experimental condition, which corresponds to the same-numbered h5ad file in the `predictions` output and the same-numbered row in the `experiments.csv` output, is used as a baseline for computing performance improvement over baseline.
- `network_datasets` describes a GRN using the same names as our network collection. The behavior is complicated because the network collection separates out tissue-specific subnetworks. The value is a dict where keys are network sources and values are (sub-)dicts controlling specific behaviors.
    - To use certain subnetworks, set `subnets` to a list naming them. To use all, set subnets to all (default).
    - To take the union of the subnetworks, set `do_aggregate_subnets` to `true`. To keep subnetworks separate, set `do_aggregate_subnets` to `false` (default).
    - You can use `empty` or `dense` for a network with no edges or all possible edges. Default is `dense`. 
    - An example from experiment `1.3.2_9`:

            "network_datasets": {
                "empty": {},
                "dense": {},
                "magnum_compendium_394": {
                    "subnets": [
                        "retinal_pigment_epithelial_cells.parquet",
                        "chronic_myelogenous_leukemia_cml_cell_line.parquet",
                        "teratocarcinoma_cell_line.parquet",
                        "lung_adenocarcinoma_cell_line.parquet",
                        "breast_carcinoma_cell_line.parquet",
                        "embryonic_kidney_cell_line.parquet",
                        "hepatocellular_carcinoma_cell_line.parquet",
                        "epitheloid_cancer_cell_line.parquet",
                        "acute_myeloid_leukemia_fab_m5_cell_line.parquet"
                    ],
                    "do_aggregate_subnets": false
                }
            }
        

There are many other keys describing the data splitting, the x/color/facet of the automated plots, the way network structures are used/pruned/ignored, and the regression methods used. Use `perturbation_benchmarking_package.experimenter.get_default_metadata()` to see the default values of each metadata field. Use `perturbation_benchmarking_package.experimenter.get_required_keys()` to see which keys are required.  Use `perturbation_benchmarking_package.experimenter.get_optional_keys()` to learn about optional keys.  
