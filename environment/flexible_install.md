
#### Flexible install

In case of different operating systems, exact conda environment reproduction is infeasible. However, you may be able to carry out some of our experiments even if some of these dependencies cannot be installed. You can try the code given below. This code requires mamba and conda to be installed already.

```bash
mamba env create --name ggrn --file environment/conda_inputs.yaml
conda activate ggrn
pip install vl-convert-python
# Avoid interfering with gimmemotifs 0.17 or with the PyTorch install by using --no-deps
# The deps should be taken care of by the above.
mamba install 'pandas>=1.2' --no-deps
pip install cell-gears==0.0.4  --no-deps
pip install celloracle==0.12.0 --no-deps
pip install prescient==0.1.0   --no-deps 
pip install geomloss==0.2.3    --no-deps 
pip install git+https://github.com/bowang-lab/scFormer@2df344a --no-deps
pip install 'scib>=1.0.3' --no-deps
for p in load_networks load_perturbations ggrn_backend2 ggrn perturbation_benchmarking_package geneformer_embeddings
do
    pip install "${p} @ git+https://github.com/ekernf01/${p}"
done
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
pip install .
cd ..


# Get data collections from Zenodo 
# accessory data, e.g. pLI and list of TF names
wget https://zenodo.org/api/files/077c17cc-c54b-4f5c-baf0-f7aca944d523/accessory_data.zip?versionId=62235c35-cdd2-4be0-b397-069d7382a32f  \
  && unzip accessory_data.zip
# perturbations 
wget https://zenodo.org/api/files/077c17cc-c54b-4f5c-baf0-f7aca944d523/perturbation_data.zip?versionId=fa60f8f8-abf7-4d5b-81e0-92ab3c6c2340 \
  && unzip perturbation_data.zip 
# networks
wget https://zenodo.org/api/files/077c17cc-c54b-4f5c-baf0-f7aca944d523/network_collection.zip?versionId=aab4e174-e575-4653-96f8-ab2273ea896c \
  && unzip network_collection.zip
```
