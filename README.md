### GGRN: benchmarking across a Grammmar of Gene Regulatory Network inference methods 

This repo contains (the beginnings of) tools for benchmarking various strategies for predicting in detail the outcome of perturbation experiments. For more context, see the [project summary](https://docs.google.com/document/d/1vvZi5c0nU3VTvKiWpEka8AtDORxJ3Ymv0ZzoFQwoDoI/edit).

### Installation

The project is written in Python with the main dependencies being gseapy, DCD-FG, scanpy + AnnData, DuckDB, PyTorch, and the usuals (joblib, numpy, pandas, scipy, etc). The environment can be reproduced via the command below. This requires a linux64 system (we use Ubuntu 20.04 or other Debian linux). 

`conda create --name ggrn --file spec_file.txt`

N.B. `conda create` does not reproduce enviroments with pip-installed packages. We avoid using any pip-installed packages; future development should do the same or revise plans for dependency management.

### Getting started

Our experiments rely on modules in the `src` folder for making predictions and evaluating them, plus [adjacent repos](https://github.com/ekernf01/perturbation_writing) for quick loading of network structures and perturbation data. Key modules are unit tested using the `unittest` framework (`src/test_<foo>.py` tests `src/foo.py`). Current development (2022 12 13) focuses on `ggrn.py`. To get started developing, pull the latest commits from the [benchmarking repo](https://github.com/ekernf01/perturbation_benchmarking). Run the unit tests in `src/test_ggrn.py`. Run some basic experiments (suggested: `test` and `1.0_1`).

### Experiments

The benchmarks in this project are composed of small, structured folders ("Experiments"). Code and metadata supporting each Experiment resides in a folder in `experiments`, and experiments can be run via `run_experiments.sh` -- but don't do this casually because it can easily take an overnight to get through everything.

- `metadata.json`: Each metadata file must have the exact same keys as one of the templates, and the values must have the same type. The most important stuff in the metadata is:
    - An explanatory "README" field
    - the associated question (from the project's [guiding questions](https://docs.google.com/document/d/1vvZi5c0nU3VTvKiWpEka8AtDORxJ3Ymv0ZzoFQwoDoI/edit#heading=h.3lbpjmchifq2))
    - ONE associated perturbation dataset (from the perturbation dataset repo)
    - ONE OR MORE associated networks (from the network collection repo). 
- `this_experiment.py`: script that runs benchmarks & makes plots unique to this experiment. Checked into git.
- `outputs`: bigger files saved as output: predicted gene expression, per-gene summaries of prediction accuracy & trainset accuracy, pickled models, and plots. Not checked into git.

### Narrative

It can be hard to understand how all these experiments relate to one another, software-wise or science-wise. Run `src/gather_all_metadata.py` for a quick summary table describing all experiments. Look on google drive for the project summary describing the guiding questions that inspired each experiment. 
