### Installation

We use Conda + Mamba + Docker to manage most dependencies. We offer either a minimal install, which may work cross-platform but lacks access to many GRN methods, or an exact install, which can reproduce all our results but only works on Ubuntu 20.04. 

### Hardware

Certain models nominally require GPU's, but we have been able to run most experiments using a CPU, sometimes by making minimal changes to Pytorch code. See the [GGRN repo](https://github.com/ekernf01/ggrn) for details on GPU requirements of specific methods. To install with GPU available, we recommend you use the exact or minimal install above; activate the environment; and then install a gpu version of PyTorch `2.x.x`.

50GB of disk space and 64GB of RAM is enough resources to run most experiments. Certain tree-based models or large datasets (Norman especially) may require more RAM. The more benchmarks you run, the more predictions are saved and the more disk space is occupied. To re-run all experiments, we would recommend 250GB disk space. 

### Minimal install

In case of different operating systems and environments, exactly reproducing results is infeasible. However, you should be able to carry out many of our experiments, or your own new experiments, even without all dependencies. Use the commands in `install_minimal.sh`. They will install python packages in a new conda environment, and they will download about 20Gb of data from Zenodo.

Some notes:

- Data download and unzip will still be running in the background after the installer finishes. It is a ~20GB download. This means the experiments may not work immediately. If you see no `network_collection` or `perturbation_data` folders, then you need to wait for the download+unzip to finish.
- We require wget, git, and conda to be installed already. If the data download doesn't work with `wget`, you can easily rephrase it to use `curl` instead, or you can download the data manually using a web browser.
- This doesn't try to install all dependencies, so some backends may be unavailable.
- This doesn't try to install docker or singularity. If you want access to containerized methods via ggrn, you need to install Docker yourself.
- This install code is written in bash and tested on Ubuntu Linux. On a Mac, the default shell is zsh, and you may need to run a bash shell (just type `bash`) for everything to work. We are not able to support Windows.

### Exact install

To reproduce our computing environment exactly, you can start with a bash shell on a clean linux box (we have tested Rocky Linux release 8.8, Ubuntu 20.04, and Ubuntu 22.04). Use these commands. This will install mamba and many python packages, and they will download about 20Gb of data from Zenodo. The one missing piece is that this doesn't try to install docker or singularity. If you want access to containerized methods via ggrn, you need to install Docker yourself.

```bash
git clone https://github.com/ekernf01/perturbation_benchmarking
source perturbation_benchmarking/environment/install.sh
```

Now, install and configure docker. This is optional; you can run many but not all experiments without it. I cannot support you in this step; if it doesn't work, you will need to go to the official Docker instructions or another source. But this worked for me on an Amazon EC2 running Ubuntu 22.04.

```bash
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
# Now configure it so you don't need sudo for every docker command.  
# https://askubuntu.com/questions/477551/how-can-i-use-docker-without-sudo
sudo usermod -aG docker $USER
sg docker -c "bash" 
```

### How to check the installation

**Warning**: data download and unzip will still be running in the background after the installer finishes. It is a ~20GB download. This means **the experiments may not work immediately.** If you see no `network_collection` or `perturbation_data` folders, then you need to wait for the download+unzip to finish.

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
pereggrn -h # see the help page
pereggrn --experiment_name "1.0_0" --amount_to_do models --no_skip_bad_runs # Run a simple benchmark
```

**Warning**: data download and unzip will still be running in the background after the installer finishes. It is a ~20GB download. This means **the experiments may not work immediately.** If you see no `network_collection` or `perturbation_data` folders, then you need to wait for the download+unzip to finish.