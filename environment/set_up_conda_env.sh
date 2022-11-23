conda create -n cell_type_grn_transfer python=3.9
conda activate  cell_type_grn_transfer
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install numba cython pybedtools jupyter jupyter-lab pysam fa2 
conda install -c conda-forge scanpy python-igraph leidenalg
conda install -c conda-forge altair vega_datasets altair_saver
conda install pyarrow
conda install -c conda-forge python-duckdb
conda install -n cell_type_grn_transfer ipykernel --update-deps --force-reinstall
conda install regex