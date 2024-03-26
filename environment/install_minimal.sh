mkdir expression_forecasting_benchmarks
cd expression_forecasting_benchmarks
# Get data collections from Zenodo 
# accessory data, e.g. pLI and list of TF names
curl https://zenodo.org/record/10436339/files/accessory_data.zip -O -s  && unzip accessory_data.zip > accessory_data.log &
# perturbations 
curl https://zenodo.org/record/10436339/files/perturbation_data.zip -O -s perturbation_data.log && unzip perturbation_data.zip > perturbation_data.log &
# networks
curl https://zenodo.org/record/10436339/files/network_collection.zip -O -s network_collection.log  && unzip network_collection.zip > network_collection.log &

# Get Python packages
mamba env create --name ggrn --file perturbation_benchmarking/environment/conda_inputs_minimal.yaml
conda activate ggrn
pip install vl-convert-python
for p in pereggrn_networks pereggrn_perturbations pereggrn ggrn ggrn_backend2 geneformer_embeddings
do
    git clone "https://github.com/ekernf01/${p}"
    pip install -e "${p}" --no-deps 
done