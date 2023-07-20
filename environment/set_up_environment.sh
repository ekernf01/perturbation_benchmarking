# This script sets up a new box (for us, usually an AWS EC2 instance) to run benchmarking analyses. 
 
# Get mamba
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh -b
source "${HOME}/mambaforge/etc/profile.d/conda.sh"

# Set up Conda env
git clone http://github.com/ekernf01/perturbation_benchmarking #currently private
cd perturbation_benchmarking
# If you have a GPU, you can use conda_list_explicit_gpu.txt.
mamba create --name ggrn --file environment/conda_list_explicit.txt
conda activate ggrn
pip install vl-convert-python
# Avoid interfering with gimmemotifs 0.17 or with the PyTorch install by using --no-deps
# The deps should be taken care of by the above.
mamba install 'pandas>=1.2' --no-deps
pip install cell-gears==0.0.4  --no-deps
pip install celloracle==0.12.0 --no-deps
pip install prescient==0.1.0   --no-deps 
pip install geomloss==0.2.3    --no-deps 
pip install git+https://huggingface.co/ctheodoris/Geneformer@f0b6641f --no-deps
pip install git+https://github.com/bowang-lab/scFormer@2df344a --no-deps
pip install 'scib>=1.0.3' --no-deps

# Install our packages
cd
for p in load_networks load_perturbations ggrn_backend2 ggrn perturbation_benchmarking_package 
do
    git clone http://github.com/ekernf01/${p}
    pip install -e $p
done

# Get data collections
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install

# For now, get data from AWS S3
# Eventually this should be a Zenodo download or something public.
aws configure 
for d in network_collection perturbation_data accessory_data  
do
    mkdir $d
    aws s3 sync s3://cahanlab/eric.kernfeld/eric_laptop/research/projects/perturbation_prediction/cell_type_knowledge_transfer/$d/ $d
done
