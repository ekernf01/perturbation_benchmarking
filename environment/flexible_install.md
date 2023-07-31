
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
for p in load_networks load_perturbations ggrn_backend2 ggrn perturbation_benchmarking_package
do
    pip install "${p} @ git+https://github.com/ekernf01/${p}"
done
git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
pip install .
cd ..
```