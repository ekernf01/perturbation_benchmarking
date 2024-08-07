mkdir expression_forecasting_benchmarks
cd expression_forecasting_benchmarks
# Get data collections from Zenodo 
# accessory data, e.g. pLI and list of TF names
wget https://zenodo.org/record/13255724/files/accessory_data.zip  && unzip accessory_data.zip &
# perturbations 
wget https://zenodo.org/record/13255724/files/perturbation_data.zip && unzip perturbation_data.zip &
# networks
wget https://zenodo.org/record/13255724/files/network_collection.zip && unzip network_collection.zip &

# Get experiment metadata and project folder layout
git clone https://github.com/ekernf01/perturbation_benchmarking
# Install python packages
conda create -n ggrn_minimal
conda activate ggrn_minimal
conda install -y pip
pip install vl-convert-python
pip install ray[tune]
for p in pereggrn_networks pereggrn_perturbations ggrn ggrn_backend2 geneformer_embeddings pereggrn
do
    git clone "https://github.com/ekernf01/${p}"
    pip install -e "${p}"
done
echo "Installation has finished, but data downloads may continue in the background, so it may not work right away."
echo "Test your installation:"
echo "    conda activate ggrn_minimal"
echo "    cd perturbation_benchmarking"
echo "    pereggrn -h # see the help page"
echo "    pereggrn --input perturbation_benchmarking --output example_output --experiment_name '1.0_0' --networks network_collection --data perturbation_data --amount_to_do models --no_skip_bad_runs"