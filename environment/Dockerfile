# from https://pythonspeed.com/articles/activate-conda-dockerfile/
FROM mambaorg/micromamba:1.4-focal
COPY --chown=$MAMBA_USER:$MAMBA_USER conda_inputs.yaml .
RUN micromamba create -y -f conda_inputs.yaml && \
    micromamba clean --all --yes

# Install python packages
SHELL ["conda", "run", "-n", "ggrn", "/bin/bash", "-c"]
RUN pip install vl-convert-python
RUN pip install cell-gears==0.0.4  --no-deps
RUN pip install celloracle==0.12.0 --no-deps
RUN pip install prescient==0.1.0   --no-deps 
RUN pip install geomloss==0.2.3    --no-deps 
RUN pip install git+https://huggingface.co/ctheodoris/Geneformer@50e921d --no-deps
RUN pip install git+https://github.com/bowang-lab/scFormer@2df344a --no-deps
RUN pip install 'scib>=1.0.3' --no-deps
RUN pip install biomart --no-deps

# Install our packages
RUN for p in load_networks load_perturbations ggrn_backend2 ggrn perturbation_benchmarking_package geneformer_embeddings \
do \
    git clone http://github.com/ekernf01/${p} && pip install -e $p --no-deps \
done

# Get accessory data, e.g. pLI and list of TF names
RUN apt-get -qq update
RUN apt-get -yq install wget
RUN wget https://zenodo.org/record/8071809/files/accessory_data.zip && unzip accessory_data.zip &
