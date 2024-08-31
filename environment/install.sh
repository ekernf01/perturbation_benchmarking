# This script sets up a new box (for us, usually an AWS EC2 instance) to run benchmarking analyses. 
 
# Get data collections from Zenodo 
sudo apt install unzip
# accessory data, e.g. pLI and list of TF names
wget https://zenodo.org/record/13345104/files/accessory_data.zip && unzip accessory_data.zip &
# perturbations 
wget https://zenodo.org/record/13345104/files/perturbation_data.zip && unzip perturbation_data.zip &
# networks
wget https://zenodo.org/record/13345104/files/network_collection.zip && unzip network_collection.zip &

# Get mamba
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh -b
source "${HOME}/mambaforge/etc/profile.d/conda.sh"

# Set up Conda env
# If you have a GPU, you can use conda_list_explicit_gpu.txt.
mamba create --name ggrn --file perturbation_benchmarking/environment/conda_list_explicit.txt
conda activate ggrn
# Why --no-deps? Makes sure every version is pinned explicitly and is compatible with the other packages.
pip install vl-convert-python==1.4.0 --no-deps
pip install git+https://github.com/snap-stanford/GEARS@df09d7a --no-deps
# PRESCIENT and CO are now used thru docker, but I am leaving this alone for backwards compatibiliy. 
pip install celloracle==0.12.0 --no-deps
pip install prescient==0.1.0   --no-deps 
pip install geomloss==0.2.3    --no-deps 
pip install git+https://github.com/bowang-lab/scFormer@2df344a --no-deps
pip install 'scib>=1.0.3' --no-deps
pip install biomart==0.9.2 --no-deps
pip install msgpack==1.0.8 --no-deps
pip install tensorboardX>=1.9 --no-deps
pip install ray[tune]==2.6.2 --no-deps
pip install scrublet==0.2.3 --no-deps
pip install pot==0.9.3 --no-deps
pip install wot==1.0.8.post2 --no-deps

# We need a specific version of Geneformer. We use `git lfs pull` because we need certain model files locally. 
echo "Cloning geneformer -- this could take a long time."
git lfs install
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
git checkout 50e921d
pip install .
git lfs pull

# Install our packages
for p in pereggrn_networks pereggrn_perturbations pereggrn ggrn ggrn_backend2 ggrn_backend3 geneformer_embeddings 
do
    git clone http://github.com/ekernf01/${p} 
    pip install -e $p --no-deps 
done

echo "The package installation has finished, but the data download and unzip may still be running in the background, so it may not work right away."
echo "Test your installation:"
echo "    conda activate ggrn"
echo "    pereggrn -h # see the help page"
echo "    pereggrn --experiment_name '1.0_0' --amount_to_do models --no_skip_bad_runs # Run a simple benchmark "
