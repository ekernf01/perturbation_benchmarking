### GGRN: benchmarking across a Grammmar of Gene Regulatory Network inference methods 

This repo contains tools for benchmarking various strategies for predicting in detail the outcome of perturbation experiments. For more context, see the [project summary](https://docs.google.com/document/d/1vvZi5c0nU3VTvKiWpEka8AtDORxJ3Ymv0ZzoFQwoDoI/edit).

### Installation

The project is written in Python with the main dependencies being gseapy, DCD-FG, scanpy + AnnData, DuckDB, PyTorch, and the usuals (joblib, numpy, pandas, scipy, etc). The environment can be exactly reproduced via the command below. This requires a linux64 system (we use Ubuntu 20.04 or other Debian linux). 

`conda create --name ggrn --file spec_file.txt`

N.B. `conda create` does not reproduce enviroments with pip-installed packages. We currently avoid using any pip-installed packages; future development should either do the same or revise plans for dependency management.

**TO DO: explain how to install our packages** Our experiments rely on [adjacent repos](https://github.com/ekernf01/perturbation_writing) for network modeling and for quick loading of network structures and perturbation data. 

### Experiments

The benchmarks in this project are composed of small, structured folders ("Experiments"). Code and metadata supporting each Experiment resides in a folder in `experiments`, and experiments can be run via `run_experiments.sh`.

- `outputs`: these results are auto-generated by `do_one_experiment.py`.
- `metadata.json`: A complete specification of the experiment, i.e. every experiment is run using the same code applied to different metadata. We do not have a formal specification for the metadata yet, but validation and defaults are done using a function from this repo's [companion package](https://github.com/ekernf01/perturbation_benchmarking_package). The most important information in the metadata is:
    - An explanatory "README" field
    - the associated question (from the project's [guiding questions](https://docs.google.com/document/d/1vvZi5c0nU3VTvKiWpEka8AtDORxJ3Ymv0ZzoFQwoDoI/edit#heading=h.3lbpjmchifq2))
    - ONE associated perturbation dataset (from the perturbation dataset repo)
    - ONE OR MORE associated networks (from the network collection repo). 

### Narrative

It can be hard to understand how all these experiments relate to one another, software-wise or science-wise. Run `gather_experiment_metadata.py` to produce a quick summary table describing all experiments. 

**TO DO: make a table specifying the experiment that each figure originated from**