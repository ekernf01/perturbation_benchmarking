# Differentiable Causal Discovery with Factor Graphs

This code is an interface to software for ["Large-Scale Differentiable Causal Discovery of Factor Graphs"](https://arxiv.org/abs/2206.07824), lightly modified from the [original repo](https://github.com/Genentech/dcdfg). That repository was originally forked from [DCDI](https://github.com/slachapelle/dcdi). Please refer to the NOTICE file for licensing information.

## Usage

    import dcdfg_wrapper
    factor_graph_model = DCDFGWrapper()
    factor_graph_model.train(my_anndata)
    predictions_anndata = factor_graph_model.predict([('POU5F1', 5), ('NANOG', 0)])
    
## Requirements

Python 3.9+ is required. To install:

    pip install -r  requirements.txt


wandb is required for now. Follow the steps [here](https://docs.wandb.ai/quickstart).