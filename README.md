### Benchmarking predictions of transcription following genetic perturbations

This repo contains (the beginnings of) tools for benchmarking various strategies for predicting in detail the outcome of perturbation experiments. For more context, see the [project summary](https://docs.google.com/document/d/1vvZi5c0nU3VTvKiWpEka8AtDORxJ3Ymv0ZzoFQwoDoI/edit).

### Environment

The project is so far written in Python with the main dependencies being Scanpy and CellOracle. Dependencies are managed by Conda (environment is called `cell_type_grn_transfer`). All development is currently done on Ubuntu 20.04, and currently, nothing is containerized. See `environment` for more details.

### Experiments

The benchmarks in this project are composed of small, structured folders ("Experiments"). Code and metadata supporting each Experiment resides in a folder in `experiments`, and experiments can be run via `run_experiments.sh` -- but don't do this casually because it can easily take an overnight to get through everything.

- `metadata.json`: experimental metadata in a strictly structured format. Each metadata file must have the exact same fields as one of the templates, and those fields must have the same types. The most important stuff in the metadata is:
    - An explanatory "README" field
    - the associated question (from the "guiding questions")
    - ONE associated perturbation dataset (from the perturbation dataset repo)
    - ONE OR MORE associated networks (from the network collection repo). 
- `this_experiment.py`: script that runs benchmarks & makes plots unique to this experiment. Checked into git.
- `outputs`: bigger files saved as output: predicted gene expression, per-gene summaries of prediction accuracy, celloracle objects, plots. Not checked into git.

### Narrative

It can be hard to understand how all these experiments relate to one another, software-wise or science-wise. Run `src/gather_all_metadata.py` for a quick summary table describing all experiments. Look on google drive for the project summary describing the guiding questions that inspired each experiment. 

### Local codebase

Our experiments rely on modules in the `src` folder for making predictions and evaluating them, plus [adjacent repos](https://github.com/ekernf01/perturbation_writing) for quick loading of network structures and perturbation data. Key modules are unit tested using the `unittest` framework (`src/test_<foo>.py` tests `src/foo.py`). Current development (2022 10 27) focuses on `predict.py`, with objectives described in `src/README.md`. 
