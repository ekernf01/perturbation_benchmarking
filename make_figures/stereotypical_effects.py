import numpy as np
import pandas as pd
import pereggrn_perturbations
import altair as alt
import re
pereggrn_perturbations.set_data_path('../../perturbation_data/perturbations')
DATASET_ORDER = [
  "nakatake",
  "nakatake\nscrna\nsimulated",
  "joung",
  "norman",
  "replogle1",
  "replogle3",
  "replogle4",
  "adamson",
  "replogle2",
  "replogle2 large effect",
  "replogle2 tf only",
  "replogle2_large_effect",
  "replogle2_tf_only",
  "replogle2\nlarge effect",
  "replogle2\nlarge\neffect",
  "replogle2\ntf only",
  "freimer",
  "dixit", 
  "frangieh_IFNg_v1",
  "frangieh\nIFNg v1",
  "frangieh\nIFNg\nv1",
  "frangieh IFNg v1",
  "frangieh_IFNg_v2",
  "frangieh\nIFNg v2",
  "frangieh\nIFNg\nv2",
  "frangieh IFNg v2",
  "frangieh_IFNg_v3",
  "frangieh\nIFNg v3",
  "frangieh\nIFNg\nv3",
  "frangieh IFNg v3"
]
datasets_used = [d for d in DATASET_ORDER if d in [
    "nakatake", 
    "replogle1",
    "replogle2",
    "replogle2\ntf only",
    "replogle2\nlarge effect",
    "replogle3",
    "replogle4",
    "joung",
    "freimer",
]]

num_pairs = 100
correlations = dict()
for dataset in datasets_used:    
    adata = pereggrn_perturbations.load_perturbation(re.sub("\n| ", "_", dataset))
    correlations[dataset] = pd.DataFrame(index = range(num_pairs), columns = ["correlation", "dataset"])
    i=0
    while i < num_pairs:
        control1 = np.random.choice(adata.obs_names[adata.obs["is_control"]])
        control2 = np.random.choice(adata.obs_names[adata.obs["is_control"]])
        treatment1 = np.random.choice(adata.obs.loc[~adata.obs["is_control"], "perturbation"])
        treatment2 = np.random.choice(adata.obs.loc[~adata.obs["is_control"], "perturbation"])
        if treatment1==treatment2:
            continue
        if control1==control2:
            continue  
        lfc1 = adata[control1].X.mean(axis=0) - adata[adata.obs["perturbation"]==treatment1].X.mean(axis=0)
        lfc2 = adata[control2].X.mean(axis=0) - adata[adata.obs["perturbation"]==treatment2].X.mean(axis=0)
        correlations[dataset].loc[i, "correlation"] = np.corrcoef(lfc1, lfc2)[0, 1]
        correlations[dataset].loc[i, "dataset"] = dataset
        i = i + 1
correlations = pd.concat(correlations.values())
box_plot = alt.Chart(correlations).mark_boxplot(color = "black").encode(
    x=alt.X('dataset:N', title='', sort = datasets_used),
    y=alt.Y('correlation:Q', title='Pearson correlation between log FC within 100 independent control-treatment pairs'),
).properties(
    title='Correlations by dataset', 
)
red_line = alt.Chart(pd.DataFrame({'y': [0]})).mark_rule(color='red').encode(
    y='y:Q'
)

# Combine the box plot and the red line
final_chart = box_plot + red_line

# Save the chart
final_chart.save("plots/stereotypical_responses.svg")