# This script sets up a new box (for us, usually an AWS EC2 instance) to run benchmarking analyses. 
 
# Get mamba
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh -b
bash #restart shell for mamba to be available

# Set up Conda env
git clone http://github.com/ekernf01/perturbation_benchmarking #currently private
cd perturbation_benchmarking
mamba create --name ggrn --file spec_file.yaml
conda activate ggrn
pip install cell-gears==0.0.2  --no-deps
pip install celloracle==0.12.0 --no-deps
pip install prescient==0.1.0   --no-deps 
pip install geomloss==0.2.3    --no-deps 
pip install vl-convert-python

# Install our packages
cd
for p in load_networks load_perturbations ggrn_backend2 ggrn_backend3 ggrn perturbation_benchmarking_package 
do
    git clone http://github.com/ekernf01/${p}
    pip install -e $p
done

# Get data collections
curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
unzip awscliv2.zip
sudo ./aws/install
aws configure #then put in pass and stuff. Eventually this should be a Zenodo download or something public.
for d in network_collection perturbation_data accessory_data
do
    mkdir $d
    aws s3 sync s3://cahanlab/eric.kernfeld/eric_laptop/research/projects/perturbation_prediction/cell_type_knowledge_transfer/$d/ $d
done