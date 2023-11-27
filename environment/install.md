### Installation

##### OS + package manager

We use Conda + Mamba to manage most dependencies. We offer either a minimal install, which may work cross-platform but lacks access to many GRN methods, or an exact install, which requires can reproduce all our results but only works on Ubuntu 20.04. 

##### Hardware

Certain models nominally require GPU's, but we have been able to run most experiments using a CPU, sometimes by making minimal changes to Pytorch code. See the [GGRN repo](https://github.com/ekernf01/ggrn) for details on GPU requirements of specific methods. 50GB of disk space and 64GB of RAM is enough resources to run most experiments. Certain tree-based models require more RAM, as does GEARS on the Norman data. The more benchmarks you run, the more predictions are saved and the more disk space is occupied. To re-run all experiments, we would recommend 250GB disk space to be safe. 

##### Desired results

For each of the options below, the desired results are a Conda environment called 'ggrn' and a folder containing four sub-folders: three data collections and the benchmark experiments.

```
├── accessory_data
├── network_collection
├── perturbation_data
├── perturbation_benchmarking 
```

### Flexible install

In case of different operating systems, exact conda environment reproduction is infeasible. However, you may be able to carry out some of our experiments, or your own new experiments, even without all dependencies. Follow the minimal installation instructions below. Some notes:

- This code requires mamba to be installed already. Install mamba using the [official instructions](https://mamba.readthedocs.io/en/latest/installation.html).
- This install code uses bash. On a Mac, the default shell is zsh, and you may need to run a bash shell (just type `bash`) for everything to work. On Windows, we are only beginning to test these instructions. Some of the major differences: bash commands must be adapted if used via PowerShell, and conda configuration may work differently. It is easy to think you have activated a conda environment, but if you have not configured the shell to alias `activate` to `activate.bat`, activation will not work.
- If the data download doesn't work with `curl`, you can easily rephrase it to use `wget` instead, or you can download the data manually from a web browser.

Here are the steps for a minimal install. 

```bash
mkdir expression_forecasting_benchmarks
cd expression_forecasting_benchmarks
# Get data collections from Zenodo 
# accessory data, e.g. pLI and list of TF names
curl https://zenodo.org/record/8071809/files/accessory_data.zip -O -s  && unzip accessory_data.zip > accessory_data.log &
# perturbations 
curl https://zenodo.org/record/8071809/files/perturbation_data.zip -O -s perturbation_data.log && unzip perturbation_data.zip && mv perturbation_data_ perturbation_data > perturbation_data.log &
# networks
curl https://zenodo.org/record/8071809/files/network_collection.zip -O -s network_collection.log  && unzip network_collection.zip > network_collection.log &

# Get Python packages
git clone https://github.com/ekernf01/perturbation_benchmarking
mamba env create --name ggrn --file perturbation_benchmarking/environment/conda_inputs_minimal.yaml
conda activate ggrn
pip install vl-convert-python
for p in load_networks load_perturbations ggrn_backend2 ggrn perturbation_benchmarking_package geneformer_embeddings
do
    git clone "https://github.com/ekernf01/${p}"
    pip install -e "${p}" --no-deps 
done
```

You can test your installation by running this.

```bash
cd perturbation_benchmarking
conda activate ggrn
python do_one_experiment.py -h # see the help page
python do_one_experiment.py --experiment_name "1.0_0" --amount_to_do models
```

### Exact reproduction

To reproduce our results, you need a clean Ubuntu 20.04 box. We recommend using AWS. This could take a long time (weeks); we ran experiments bit by bit over a long period. 

The repo is under active development and may not be entirely stable or may not exactly reproduce our preprint. A list of commit hashes used for version one of our preprint can be found in the `environment` folder, and we plan to make code releases for future preprint versions or journal submissions.

```bash
git clone https://github.com/ekernf01/perturbation_benchmarking
source environment/install.sh
cd perturbation_benchmarking
source run_experiments.sh &  #This could take weeks!
Rscript make_figures/*.R 
```

Explanation:

- Run  `environment/install.sh` to install the Conda environment and download the data.
- Experiments can be run via `source run_experiments.sh &`. This takes a long time. You can also run individual experiments via `do_one_experiment.py`.  
- Once the experiments are done, figures can be produced using the R scripts in `make_figures`. This requires some common basic R packages like ggplot2 that are not included in our environment setup -- let us know if you have trouble installing them.


