conda create -n cell_type_grn_transfer python=3.9
conda activate  cell_type_grn_transfer
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install numba cython pybedtools jupyter jupyter-lab pysam fa2 
conda install -c conda-forge leidenalg=0.8.10
pip3 install git+https://github.com/morris-lab/CellOracle.git
conda install scanpy
conda install dask
