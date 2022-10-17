### Benchmarking predictions of transcription following genetic perturbations

This repo contains (the beginnings of) tools for benchmarking various strategies for predicting in detail the outcome of perturbation experiments. For more information, see the [project rundown](https://docs.google.com/document/d/1vvZi5c0nU3VTvKiWpEka8AtDORxJ3Ymv0ZzoFQwoDoI/edit).

### Environment

The project is so far written in Python with the main dependencies being Scanpy and CellOracle. Dependencies are managed by Conda (environment is called `cell_type_grn_transfer`). See `envrionment` for more details.

### Experiments

The benchmarks in this project are composed of small, structured folders ("Experiments"). Code and metadata supporting each Experiment resides in a folder in `experiments`.

- `metadata.json`: experimental metadata in a strictly structured format. Each metadata file must have the exact same fields as one of the templates, and those fields must have the same types. The most important stuff in the metadata is:
    - the associated question (from the "guiding questions")
    - ONE associated perturbation dataset (from the perturbation dataset repo)
    - ONE OR MORE associated networks (from the network collection repo). 
- `this_experiment.py`: script that runs benchmarks & makes plots unique to this experiment. Checked into git.
- `outputs`: bigger files saved as output: predicted gene expression, per-gene summaries of prediction accuracy, celloracle objects, plots. Not checked into git.

