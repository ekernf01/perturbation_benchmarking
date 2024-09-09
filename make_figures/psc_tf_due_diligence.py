import scanpy as sc
import pandas as pd
import numpy as np
import pereggrn_perturbations
from pereggrn import experimenter
from scipy.stats import rankdata
pereggrn_perturbations.set_data_path("../../perturbation_data/perturbations")
distance_to_targets = pd.read_csv("boyer2005targets.csv")[["GENE", "SOX2", "NANOG", "E2F4"]]
distance_to_targets = 
predictions = sc.read_h5ad("../experiments/1.0_1/outputs/predictions/7.h5ad") # condition 7 is LassoCV
observed = pereggrn_perturbations.load_perturbation("nakatake")
observed = experimenter.averageWithinPerturbation(observed)
top_genes = {}
targets = list(set(distance_to_targets["GENE"].unique()).intersection(predictions.var_names))
overlap = pd.DataFrame(index = targets, columns = ["predicted_and_observed", "predicted_and_boyer", "boyer_and_observed"])
for target in targets:
    top_genes[target] = {}
    predicted_logfc     = (predictions[                 :, target].X.mean(axis=1) - observed[observed.obs["is_control"], target].X.mean())
    observed_logfc      = (observed[predictions.obs_names, target].X.mean(axis=1) - observed[observed.obs["is_control"], target].X.mean())
    logfc = pd.DataFrame(
        {
            "observed_logfc": observed_logfc,
            "predicted_logfc": predicted_logfc,
            "observed_absolute_logfc": np.abs(observed_logfc),
            "predicted_absolute_logfc": np.abs(predicted_logfc),
        },
        index = predictions.obs_names, 
    )
    top_genes[target]["predicted"] = logfc.query("@rankdata(-predicted_absolute_logfc)<=5").index
    top_genes[target]["observed"]  = logfc.query("@rankdata(-observed_absolute_logfc)<=5").index
    top_genes[target]["boyer"] = distance_to_targets.query("GENE==@target").T[[False, True, True, True]].set_axis(['distance'], axis = 1).query("distance!='-'").index.values
    overlap.loc[target, "observed_top"] = top_genes[target]["observed"][0]
    try:
        overlap.loc[target, "predicted_top"] = top_genes[target]["predicted"][0]
    except IndexError:
        overlap.loc[target, "predicted_top"] = ""
    overlap.loc[target, "predicted"]              = len(top_genes[target]["predicted"])
    overlap.loc[target, "boyer"]                  = len(top_genes[target]["boyer"])
    overlap.loc[target, "observed"]               = len(top_genes[target]["observed"])
    overlap.loc[target, "predicted_and_observed"] = len(set(top_genes[target]["predicted"]).intersection(top_genes[target]["observed"]))
    overlap.loc[target, "predicted_and_boyer"]    = len(set(top_genes[target]["predicted"]).intersection(top_genes[target]["boyer"]))
    overlap.loc[target, "boyer_and_observed"]     = len(set(top_genes[target]["boyer"]).intersection(top_genes[target]["observed"]))


[overlap.value_counts(c) for c in overlap.columns]

