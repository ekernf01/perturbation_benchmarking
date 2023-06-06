import pandas as pd
import altair as alt
import os 

X = {}
for experiment in ["1.6.1_1", "1.6.1_3", "1.6.1_6"]:
    X[experiment] = pd.read_parquet(f"../experiments/{experiment}/outputs/evaluationPerPert.parquet")
X = pd.concat(X.values())
X.reset_index(inplace=True)

# Create the Altair chart
alt.data_transformers.disable_max_rows()
X["chart_x"] = X["regression_method"] + "_" + X["starting_expression"]
for metric in ["mae_benefit"]:
    chart = alt.Chart(X).mark_boxplot().encode(
        y=alt.Y("mae_benefit", title="MAE improvement over baseline"),
        x=alt.X('chart_x', title='Method or network'),
        column=alt.Column('perturbation_dataset', title='Perturbation dataset'),
        color=alt.Color('starting_expression', title='Starting expression'),
    ).properties(
        width=300,  # Adjust the width as needed
        height=300  # Adjust the height as needed
    ).resolve_scale(
        x='independent',
        y='independent',
    )
    os.makedirs("plots", exist_ok=True)
    chart.show()
    chart.save(f'plots/fig_dcdfg_{metric}.html', embed_options={'renderer':'svg'})

