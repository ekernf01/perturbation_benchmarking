import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import pereggrn_perturbations
pereggrn_perturbations.set_data_path('../../perturbation_data/perturbations')
import sys 
import altair as alt
import os

sys.path.append("../../perturbation_data/setup/")
import ingestion
import global_effects

def decompose_variance(adata, gene_perturbed, dataset):
    try:
        depth = np.random.choice( a = np.array(adata.raw.X.sum(1)).reshape(-1), size = 1000, replace = True )
        pre_log_scale = np.expm1(adata[0,:].X).sum()
        fraction_of_rna_mapping_to_this_gene = np.expm1(adata[:,gene_perturbed].X).mean()/pre_log_scale
        poisson_raw_counts = np.array([np.random.poisson( lam = l ) for l in fraction_of_rna_mapping_to_this_gene*depth])
        resampled = np.log1p(pre_log_scale*(poisson_raw_counts / depth))
        poisson = np.var(resampled)
    except:
        poisson = np.nan
    control = np.var(ingestion.try_toarray(adata[adata.obs["is_control"],gene_perturbed].X))
    others  = np.var(ingestion.try_toarray(adata[adata.obs["perturbation"]!=gene_perturbed,gene_perturbed].X))
    this    = np.var(ingestion.try_toarray(adata[adata.obs["is_control"] | (adata.obs["perturbation"]==gene_perturbed),gene_perturbed].X))
    return pd.DataFrame({
        "poisson": poisson, 
        "control": control, 
        "others": others, 
        "this": this,
        "gene_perturbed": gene_perturbed, 
        "dataset": dataset,
    }, index = [0])

os.makedirs("variance_decomposition", exist_ok = True)
all_variance_decomposition = []
for dataset in [
    'nakatake',
    'freimer',
    'replogle',
    'replogle2',
    'replogle3',
    'replogle4',
    'frangieh_IFNg_v1',
    'frangieh_IFNg_v2',
    'frangieh_IFNg_v3',
    'dixit',
    'adamson',
    'norman',
]:
    print(dataset)
    try:
        variance_decomposition = pd.read_csv(f"variance_decomposition/{dataset}.csv")
    except FileNotFoundError:
        adata = pereggrn_perturbations.load_perturbation(dataset)
        pt = adata.obs["perturbation_type"][0]
        variance_decomposition = []
        for gene_perturbed in adata.uns["perturbed_and_measured_genes"]:
            variance_decomposition.append(decompose_variance(adata, gene_perturbed, dataset))
        variance_decomposition = pd.concat(variance_decomposition)
        variance_decomposition.to_csv(f"variance_decomposition/{dataset}.csv")
    all_variance_decomposition.append(variance_decomposition)
all_variance_decomposition = pd.concat(all_variance_decomposition)
all_variance_decomposition

all_variance_decomposition = pd.melt(
    all_variance_decomposition, 
    id_vars=['gene_perturbed', 'dataset'], 
    value_vars=['poisson', 'control', 'others', 'this'], 
    var_name='source_of_variance', 
    value_name='variance'
)

chart = alt.Chart(all_variance_decomposition.groupby(["dataset", "source_of_variance"])["variance"].mean()).mark_point().encode(
    x='dataset:N',  # :N denotes a nominal (discrete) variable
    y='variance:Q',  # :Q denotes a quantitative (continuous) variable
    color='source_of_variance:N'  # Color by source_of_variance, also discrete
).properties(
    title='Scatter plot of Variance by Dataset and Source of Variance'
)

chart.display()