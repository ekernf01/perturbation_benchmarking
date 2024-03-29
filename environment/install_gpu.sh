# This script sets up a new box (for us, usually an AWS EC2 instance) to run benchmarking analyses. 
 
# Get data collections from Zenodo 
sudo apt install unzip
# accessory data, e.g. pLI and list of TF names
wget https://zenodo.org/record/10436339/files/accessory_data.zip && unzip accessory_data.zip &
# perturbations 
wget https://zenodo.org/record/10436339/files/perturbation_data.zip && unzip perturbation_data.zip &
# networks
wget https://zenodo.org/record/10436339/files/network_collection.zip && unzip network_collection.zip &

# Get mamba
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh -b
source "${HOME}/mambaforge/etc/profile.d/conda.sh"

# Set up Conda env
# If you have a GPU, you can use conda_list_explicit_gpu.txt.
mamba create --name ggrn --file perturbation_benchmarking/environment/conda_list_explicit.txt
conda activate ggrn
pip install vl-convert-python
# Avoid interfering with gimmemotifs 0.17 or with the PyTorch install by using --no-deps
# The deps should be taken care of by the above. 
pip install git+https://github.com/snap-stanford/GEARS@df09d7a --no-deps
pip install geomloss==0.2.3    --no-deps 
pip install git+https://huggingface.co/ctheodoris/Geneformer@50e921d --no-deps
pip install 'scib>=1.0.3' --no-deps
pip install biomart --no-deps
pip install ray[tune]==2.6.2
pip install scrublet --no-deps

# Install our packages
for p in pereggrn_networks pereggrn_perturbations pereggrn ggrn ggrn_backend2 geneformer_embeddings 
do
    git clone http://github.com/ekernf01/${p} 
    pip install -e $p --no-deps 
done

echo "The package installation has finished, but the data download and unzip may still be running in the background."