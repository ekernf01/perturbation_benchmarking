import pandas as pd
import altair as alt
import os 

# Figure <gears>
X = {}
for experiment in [f"1.4.2_{i}" for i in range(1, 7) if i!=4]:
    X[experiment] = pd.read_parquet(f"../experiments/{experiment}/outputs/evaluationPerPert.parquet")
X = pd.concat(X.values())
X.reset_index(inplace=True)

# Create the Altair chart
alt.data_transformers.disable_max_rows()
X["chart_x"] = X["regression_method"].astype(str) + "_using_" + X["eligible_regulators"].astype(str)
for metric in ["mae", "mse_top_20"]:
    chart = alt.Chart(X).mark_boxplot().encode(
        y=alt.Y(metric, title=metric),
        x=alt.X('chart_x', title='Method and gene-set for GO graph'),
        column=alt.Column('perturbation_dataset', title='Perturbation Dataset'),
        row=alt.Row('desired_heldout_fraction', title='Desired Heldout Fraction'),
        color=alt.Color('eligible_regulators', title='Eligible Regulators')
    ).properties(
        width=300,  # Adjust the width as needed
        height=300  # Adjust the height as needed
    )
    os.makedirs("plots", exist_ok=True)
    chart.save(f'plots/fig_gears_{metric}.html', embed_options={'renderer':'svg'})

