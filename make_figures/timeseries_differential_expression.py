import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import pereggrn_perturbations
pereggrn_perturbations.set_data_path('../../perturbation_data/perturbations')
import sys 
import altair as alt
import os

dataset = 'definitive_endoderm'
adata = pereggrn_perturbations.load_perturbation(dataset,   is_timeseries=True)
adata.obs['timepoint'] = adata.obs['timepoint'].astype('str')
adata.uns["log1p"]["base"] = 2
sc.tl.rank_genes_groups(adata, "timepoint", groups=['4'], reference='0', method='wilcoxon')
human_tf = pd.read_csv('../../accessory_data/tf_lists/human.txt', header=None)
X = sc.get.rank_genes_groups_df(adata, group='4')
X = X.query("names in @human_tf[0].values")
X = X.query("pvals_adj < 0.05")
X.sort_values('scores', ascending=False, inplace=True)
X.head(30).to_csv("timeseries_plots/top30_differential_expression_endoderm.csv")

dataset = 'fantom4'
adata = pereggrn_perturbations.load_perturbation(dataset,   is_timeseries=True)
adata.obs['timepoint'] = adata.obs['timepoint'].astype('str')
sc.tl.rank_genes_groups(adata, "timepoint", groups=['96.0'], reference='0.0', method='wilcoxon')
human_tf = pd.read_csv('../../accessory_data/tf_lists/human.txt', header=None)
X = sc.get.rank_genes_groups_df(adata, group='96.0')
X = X.query("names in @human_tf[0].values")
X.sort_values('logfoldchanges', ascending=False, inplace=True)
X.head(30).to_csv("timeseries_plots/top30_differential_expression_fantom4.csv")

dataset = 'paul1'
adata = pereggrn_perturbations.load_perturbation(dataset,   is_timeseries=True)
mouse_tf = pd.read_csv('../../accessory_data/tf_lists/mouse.txt', header=None)
cell_type_groups = {
    "ME": [
        "MEP",
        "Erythroids",
        "Megakaryocytes",
        ],
    "DC": [
        "DC",
    ],
    "GM": [
        "GMP",
        "late_GMP",
        "Monocytes",
        "Granulocytes",
    ]
}
for group, celltypes in cell_type_groups.items():
    sc.tl.rank_genes_groups(adata, "cell_type", groups=celltypes, reference='rest', method='wilcoxon')
    X = sc.get.rank_genes_groups_df(adata, group=celltypes)
    X = X.query("names in @mouse_tf[0].values")
    X = X.query("pvals_adj < 0.05")
    X.sort_values('logfoldchanges', ascending=False, inplace=True)
    X.head(30).to_csv(f"timeseries_plots/top30_differential_expression_paul_{group}.csv")