import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import pereggrn_perturbations
pereggrn_perturbations.set_data_path('../../perturbation_data/perturbations')
import sys 
import altair as alt
from scipy.stats import rankdata as rank

sys.path.append("../../perturbation_data/setup/") # access our data ingestion module, which is not currently pip-installable
import ingestion
import global_effects
effects = []


DATASET_ORDER = [
    "nakatake",
    "joung",
    "norman",
    "replogle1",
    "replogle3",
    "replogle4",
    "adamson",
    "replogle2",
    "freimer",
    "dixit", 
    "frangieh_IFNg_v2",
]

# Top-n genes consistency across replicates
for dataset in ["nakatake", "freimer", "frangieh_IFNg_v3", "replogle1"]:    
    print(dataset)
    adata = pereggrn_perturbations.load_perturbation(dataset)
    baseline = adata.X[adata.obs["is_control"], :].mean(axis = 0)
    intersection = {}
    union = {}
    jaccard = {}
    for n in [20, 100, 200]:
        for perturbation in adata.obs["perturbation"].unique():
            intersection[perturbation] = set(adata.var_names)
            union[perturbation] = set()
            for i in adata.obs.query("perturbation == @perturbation").index:
                if adata.obs["is_control"][i]:
                    continue
                logfc = adata[i, :].X - baseline
                top_n_genes = set(adata.var_names[rank(-np.abs(logfc)) <= n]).copy()
                intersection[perturbation] = intersection[perturbation].intersection(top_n_genes)
                union[perturbation] = union[perturbation].union(top_n_genes)
                jaccard[perturbation] = len(intersection[perturbation]) / len(union[perturbation])
        print(n)
        print(np.array([x for x in jaccard.values()]).mean())

# Effect size and direction
for dataset in DATASET_ORDER:    
    print(dataset)
    adata = pereggrn_perturbations.load_perturbation(dataset)
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
    global_effects.quantifyEffect(adata, fname = fname, withDEG = False, withMI = False, pseudocount = 1)
    obs = adata.obs[['perturbation', 'is_control', 'logFC', 'logFCNorm2', 'logFCMean', 'expression_level_after_perturbation', 'perturbation_type']].copy()
    
    obs.loc[:, "dataset"] = dataset
    effects.append(obs)

effects = pd.concat(effects)
effects = effects.query("~is_control")
effects = effects.query("logFC != -999")
effects["guide"] = 0
alt.data_transformers.disable_max_rows()
chart = alt.Chart(effects).transform_density(
    density='logFC',  
    groupby=['dataset', 'perturbation_type'], 
    as_=['logFC', 'density'],  
    extent=[effects['logFC'].min(), effects['logFC'].max()], 
    counts=False  
).mark_area(opacity=0.75).encode(
    x=alt.X('logFC:Q', title='logFC'), 
    y=alt.Y('density:Q', title='Density'),
    color=alt.Color('dataset:N', scale=alt.Scale(domain=DATASET_ORDER, scheme = "dark2"))  # replace with your desired order
).properties(
    width=200,
    height=60
)
vline = alt.Chart(effects).mark_rule(color='black').encode(
    x='guide:Q'
)
chart = (chart + vline).facet(
    row="perturbation_type:N",
)
chart.save('plots/fig_effects.svg')


chart = alt.Chart(effects).mark_circle(size=10).encode(
    x=alt.X('logFC:Q', title='logFC of perturbed gene\'s RNA'), 
    y=alt.Y('logFCMean:Q', title='Mean absolute logFC, all genes'),
    color=alt.Color('dataset:N', scale=alt.Scale(domain=DATASET_ORDER, scheme = "dark2"))  # replace with your desired order
).properties(
    width=200,
    height=100
)
chart = chart.facet(
    row="perturbation_type:N",
).resolve_scale(y='independent') 
chart.save('plots/fig_effects2.svg')

