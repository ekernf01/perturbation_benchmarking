import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import load_perturbations
load_perturbations.set_data_path('../../perturbation_data/perturbations')
import sys 
import altair as alt

sys.path.append("../../perturbation_data/setup/")
import ingestion
import global_effects
effects = []
for dataset in [
    'nakatake',
    'freimer',
    'replogle',
    'replogle2',
    'replogle3',
    'replogle4',
    'frangieh_IFNg_v2',
    'dixit',
    'adamson',
    'norman',
]:    
    print(dataset)
    adata = load_perturbations.load_perturbation(dataset)
    pt = adata.obs["perturbation_type"][0]
    uns = adata.uns.copy()
    try:
        if adata.X.sum() == adata.raw.X.sum(): # We filled in log1p normalized data into the .raw slot for datasets obtained from GEARS. 
            adata.raw = anndata.AnnData(X = np.exp(adata.raw.X.toarray()) - 1)
        adata = ingestion.aggregate_by_perturbation(adata, group_by = ["perturbation"], use_raw = True)
        sc.pp.normalize_total(adata)
    except:
        pass
    adata.uns = uns
    adata = ingestion.describe_perturbation_effect(adata, perturbation_type = pt)
    consistency = ingestion.checkConsistency(adata, pt)
    adata.obs["logFC"] = consistency[1]
    print("Consistency:")
    print(pd.Series(consistency[0]).value_counts()) 
    fname = "global_effects/" + dataset + ".txt"
    global_effects.quantifyEffect(adata, fname = fname, withDEG = False)
    obs = adata.obs[['perturbation', 'is_control', 'logFC', 'logFCNorm2', 'logFCMean', 'expression_level_after_perturbation', 'perturbation_type']].copy()
    
    obs.loc[:, "dataset"] = dataset
    effects.append(obs)

effects = pd.concat(effects)
effects = effects.query("~is_control")
effects = effects.query("logFC != -999")


chart = alt.Chart(effects).transform_density(
    density='logFC',  
    groupby=['dataset'], 
    as_=['logFC', 'density'],  
    extent=[effects['logFC'].min(), effects['logFC'].max()], 
    counts=False  
).mark_area(opacity=0.5).encode(
    x=alt.X('logFC:Q', title='logFC'), 
    y=alt.Y('density:Q', title='Density'),
    color='dataset:N'
).properties(
    width=280,
    height=220
)
alt.data_transformers.disable_max_rows()
chart.save('plots/fig_effects.svg')


chart = alt.Chart(effects).mark_circle(size=10).encode(
    x=alt.X('logFC:Q', title='logFC of perturbed gene\'s RNA'), 
    y=alt.Y('logFCMean:Q', title='Mean absolute logFC, all genes'),
    color=alt.Color('dataset:N') 
).properties(
    width=280,
    height=220
)
chart.save('plots/fig_effects2.svg')
