 ## A systematic comparison of computational methods for expression forecasting with [PEREGGRN](https://github.com/ekernf01/pereggrn) 
 
This repo contains benchmark experiments to evaluate various strategies for predicting gene expression after knockout, knockdown, or overexpression. 

![image](https://github.com/ekernf01/perturbation_benchmarking/assets/5271803/ae7a5c86-dca6-49be-b048-743f8e110a18)

- For context and key results, see our [preprint](https://www.biorxiv.org/content/10.1101/2023.07.28.551039v2). 
- To interact with the evaluation results and see the source data for our figures, [download (300MB)](https://zenodo.org/records/13785280/files/evaluation_results.zip?download=1) them from zenodo (DOI: 10.5281/zenodo.13785280). 
- To repeat our experiments or run your own, see the [pereggrn](https://github.com/ekernf01/pereggrn) benchmarking software ([tutorial](https://github.com/ekernf01/pereggrn/blob/main/docs/tutorial.md), how to [add your own method](https://github.com/ekernf01/pereggrn/blob/main/docs/how_to.md#how-to-evaluate-a-new-method)).
- If there's something you cannot find, go ahead and file a github issue -- with your input, we hope to improve the project.

### Related infrastructure

This project is tightly coupled with our collections of data, our GGRN package for dynamic models of gene regulatory networks, and our PEREGGRN package containing benchmarking infrastructure. More information:

- Perturbation data, the network collection, and some accessory data (e.g. a list of TF's) are on Zenodo with DOI `10.5281/zenodo.10436339`.
    - Our code expects each of those three folders to be unzipped and placed adjacent to this repo.
    - Use our [perturbation loader](https://github.com/ekernf01/pereggrn_perturbations) and [network loader](https://github.com/ekernf01/pereggrn_networks) to easily access and validate data from Python.
- [GGRN](https://github.com/ekernf01/ggrn), the Grammar of Gene Regulatory Networks, offers flexible combination of different features for regulatory network inference.
- [PEREGGRN](https://github.com/ekernf01/pereggrn), PErturbation Response Evaluation via a Grammar of Gene Regulatory Networks, helps conduct the experiments that are specified in this repo.
- Certain additional experiments are implemented in [our fork of DCD-FG](https://github.com/ekernf01/dcdfg).


