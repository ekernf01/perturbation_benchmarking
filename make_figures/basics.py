import pandas as pd
import altair as alt
import os 

X = {}
for experiment in ["1.0_2", "1.0_3", "1.4.5_1", "1.4.7_1"]:
    X[experiment] = pd.read_parquet(f"../experiments/{experiment}/outputs/evaluationPerPert.parquet")
X = pd.concat(X.values())
X.reset_index(inplace=True)

# Create the Altair chart
alt.data_transformers.disable_max_rows()
X["chart_x"] = [X.loc[a,b] for a,b in zip(X.index, X["factor_varied"])]
for metric in ["mae_benefit"]:
    chart = alt.Chart(X).mark_boxplot().encode(
        y=alt.Y("mae_benefit", title="MAE improvement over baseline"),
        x=alt.X('chart_x', title='Method or network'),
        column=alt.Column('perturbation_dataset', title='Perturbation dataset'),
        row=alt.Row('factor_varied', title='Factor varied'),
    ).properties(
        width=300,  # Adjust the width as needed
        height=300  # Adjust the height as needed
    ).resolve_scale(
        x='independent',
        y='independent',
    )
    os.makedirs("plots", exist_ok=True)
    chart.save(f'plots/fig_basics_{metric}.html', embed_options={'renderer':'svg'})

