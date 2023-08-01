
#### Flexible install

In case of different operating systems, exact conda environment reproduction is infeasible. However, you may be able to carry out some of our experiments even if some of these dependencies cannot be installed. You can try the code given below. This code requires mamba and conda to be installed already.

```bash

# Get data collections from Zenodo 
# accessory data, e.g. pLI and list of TF names
wget https://zenodo.org/record/8071809/files/accessory_data.zip && unzip accessory_data.zip &
# perturbations 
wget https://zenodo.org/record/8071809/files/perturbation_data.zip && unzip perturbation_data.zip && mv perturbation_data_ perturbation_data &
# networks
wget https://zenodo.org/record/8071809/files/network_collection.zip && unzip network_collection.zip &

# Get Python packages
mamba env create --name ggrn --file environment/conda_inputs.yaml
conda activate ggrn
pip install vl-convert-python
# Avoid interfering with gimmemotifs 0.17 or with the PyTorch install by using --no-deps
# The deps for all this should be taken care of by the above.
mamba install 'pandas>=1.2' --no-deps
pip install cell-gears==0.0.4  --no-deps
pip install celloracle==0.12.0 --no-deps
pip install prescient==0.1.0   --no-deps 
pip install geomloss==0.2.3    --no-deps 
pip install git+https://github.com/bowang-lab/scFormer@2df344a --no-deps
pip install 'scib>=1.0.3' --no-deps
for p in load_networks load_perturbations ggrn_backend2 ggrn perturbation_benchmarking_package geneformer_embeddings
do
    git clone "https://github.com/ekernf01/${p}"
    pip install -e "${p}" --no-deps 
done
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
pip install .
cd ..
```
