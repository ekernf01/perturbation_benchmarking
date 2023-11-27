## A systematic comparison of computational methods for expression forecasting

This repo contains benchmark experiments to evaluate various strategies for predicting gene expression after knockout, knockdown, or overexpression. 

![image](https://github.com/ekernf01/perturbation_benchmarking/assets/5271803/ae7a5c86-dca6-49be-b048-743f8e110a18)

- For context and key results, see our [preprint](https://www.biorxiv.org/content/10.1101/2023.07.28.551039v1).
- For a summary of our guiding questions and experiments, see `all_experiments.tsv` and `guiding_questions.txt`. 
- For the installation guide: see [`environment/install.md`](https://github.com/ekernf01/perturbation_benchmarking/blob/main/environment/install.md).
- For a reference describing all possible inputs and outputs: see [docs/reference.md](https://github.com/ekernf01/perturbation_benchmarking/blob/main/docs/reference.md).
- For how-to recipes to design your own experiments, see [docs/how_to.md](https://github.com/ekernf01/perturbation_benchmarking/blob/main/docs/how_to.md).
- If there's something you cannot find, go ahead and file a github issue -- with your input, we hope to improve the project.

### Related infrastructure

This project is tightly coupled with our collections of data, our GGRN package for network inference, and a companion package containing benchmarking infrastructure. More information:

- Perturbation data, the network collection, and some accessory data (e.g. a list of TF's) are on Zenodo with DOI `10.5281/zenodo.8071809`.
    - Our code expects each of those three folders to be unzipped and placed adjacent to this repo.
    - Use our [perturbation loader](https://github.com/ekernf01/load_perturbations) and [network loader](https://github.com/ekernf01/load_networks) to easily access and validate data from Python.
- [GGRN](https://github.com/ekernf01/ggrn) offers flexible combination of different features for regulatory network inference.
- Our [perturbation benchmarking package](https://github.com/ekernf01/perturbation_benchmarking_package) helps conduct the `Experiment`s that are specified in this repo.
- Certain additional experiments are implemented in [our fork of DCD-FG](https://github.com/ekernf01/dcdfg).


