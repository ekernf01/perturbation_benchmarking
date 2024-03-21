## March 2024 updates to the benchmarking framework for timeseries prediction

- Want to be able to make predictions for any combo of `cell_type`, `timepoint`, `perturbation`, `expression_level_after_perturbation`, `prediction_timescale`. 
- For existing experiments, default is one prediction per test sample, with different `expression_level_after_perturbation`, `perturbation` but same `cell_type`, `timepoint`, `prediction_timescale`
- For TS, default will be only one `expression_level_after_perturbation` but all `cell_type`, `timepoint`, `perturbation`, `prediction_timescale`
- The interface between GGRN and benchmarking code will have to change -- cannot currently carry all this info. 

## Desired new interface: 

#### Pereggrn call to ggrn

- in GRN.predict(), the `perturbations` arg will be replaced by a new arg called `predictions_metadata`. It will be a pandas dataframe with columns `cell_type`, `timepoint`, `perturbation`, and `expression_level_after_perturbation`. 
    - It will default to `predictions.obs` or `starting_expression.obs` if those are provided.
        - **TODO**: check compatibility of sizes; check that either is provided; warn if both are provided
    - It will be provided to Dockerized methods as a csv. 
        - **TODO**: update Dock-umentation, update train.py, update ggrn.api
    - The `perturbations` arg will be completely removed.
        - **TODO**: remove this from ggrn.api, from ggrn_backend2.api, and from do_one_experiment.py
    - The meaning is "predict expression in `cell_type` at `time_point` if `perturbation` were set to `expression_level_after_perturbation`". However, methods may ignore some columns. For all aim 2 benchmarks, columns `cell_type`, `timepoint` will be placeholders and will be ignored. For GEARS, `expression_level_after_perturbation` is ignored. 

#### GGRN call to backends

- `predictions_metadata` will be passed through, except beforehand, `prediction_timescale` will be cartesian-producted with it, yielding one more column called `predictions_metadata`. To maintain DRY, any representation of `prediction_timescale` outside `predictions_metadata` should be ignored by individual backends or even better not passed to them. 
    - **TODO**: document this ... somewhere??

#### Code fragments in transit

```python
train_data_tp_by_ct = train.obs[["timepoint", "cell_type"]].unique()
combinations = list(itertools.product(train_data_tp_by_ct, perturbations, ggrn_args["prediction_timescale"]))
predictions_metadata = pd.DataFrame(
    [[c[0][0], c[0][1], c[1][0], c[1][0], c[0][0]] for c in combinations],
    columns=['timepoint', 'cell_type', 'perturbation', "expression_level_after_perturbation", 'num_steps'],
)
```