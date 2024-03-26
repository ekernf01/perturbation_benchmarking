### Installation

We use Conda + Mamba to manage most dependencies. We offer either a minimal install, which may work cross-platform but lacks access to many GRN methods, or an exact install, which requires can reproduce all our results but only works on Ubuntu 20.04. 

### Hardware

Certain models nominally require GPU's, but we have been able to run most experiments using a CPU, sometimes by making minimal changes to Pytorch code. See the [GGRN repo](https://github.com/ekernf01/ggrn) for details on GPU requirements of specific methods. 50GB of disk space and 64GB of RAM is enough resources to run most experiments. Certain tree-based models or large datasets (Norman especially) may require more RAM. The more benchmarks you run, the more predictions are saved and the more disk space is occupied. To re-run all experiments, we would recommend 250GB disk space. 

### Minimal install

In case of different operating systems, exact conda environment reproduction is infeasible. However, you may be able to carry out some of our experiments, or your own new experiments, even without all dependencies. Use these commands. They will install software via apt and mamba, and they will download about 20Gb of data from Zenodo.

```bash
git clone https://github.com/ekernf01/perturbation_benchmarking
source perturbation_benchmarking/environment/install_minimal.sh
```

Some notes:

- If the data download doesn't work with `curl`, you can easily rephrase it to use `wget` instead, or you can download the data manually from a web browser.
- We require git and mamba to be installed already. Install mamba using the [official instructions](https://mamba.readthedocs.io/en/latest/installation.html).
- This install code is written in bash. On a Mac, the default shell is zsh, and you may need to run a bash shell (just type `bash`) for everything to work. On Windows, we are only beginning to test these instructions. Some of the major differences:
    - bash commands will not work in PowerShell. The Windows Subsystem for Linux may be a good alternative.
    - conda configuration may work differently. 

### Exact install

To reproduce our computing environment exactly, you need a bash shell on a clean Ubuntu 20.04 box. We recommend using AWS. Use these commands. They will install software via apt and mamba, and they will download about 20Gb of data from Zenodo.

```bash
git clone https://github.com/ekernf01/perturbation_benchmarking
source perturbation_benchmarking/environment/install.sh
```

### How to check the installation

**Warning**: data download and unzip will still be running in the background after the installer finishes. It is a ~20GB download. This means the experiments may not work immediately.

The installation should create a Conda environment called 'ggrn' and several folders in your working directory. At minimum, there will be three data collections and the benchmark experiments.

```
├── accessory_data
├── network_collection
├── perturbation_data
├── perturbation_benchmarking 
```

You can test your installation by running this.

```bash
cd perturbation_benchmarking
conda activate ggrn
python do_one_experiment.py -h # see the help page
python do_one_experiment.py --experiment_name "1.0_0" --amount_to_do models --no_skip_bad_runs # Run a simple benchmark
```
