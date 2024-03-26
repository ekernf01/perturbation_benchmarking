## March 2024 updates to the benchmarking framework for timeseries prediction

- Want to be able to make predictions for any combo of `cell_type`, `timepoint`, `perturbation`, `expression_level_after_perturbation`, `prediction_timescale`. 
- For existing experiments, default is one prediction per test sample, with different `expression_level_after_perturbation`, `perturbation` but same `cell_type`, `timepoint`, `prediction_timescale`
- For TS, default will be only one `expression_level_after_perturbation` but all `cell_type`, `timepoint`, `perturbation`, `prediction_timescale`
- The interface between GGRN and benchmarking code will have to change -- cannot currently carry all this info. 

## Desired new interface: 

#### Pereggrn call to ggrn

- in GRN.predict(), the `perturbations` arg will be replaced by a new arg called `predictions_metadata`. It will be a pandas dataframe with columns `cell_type`, `timepoint`, `perturbation`, and `expression_level_after_perturbation`. 
    - It will default to `predictions.obs` or `starting_expression.obs` if those are provided.
    - It will be provided to Dockerized methods as a csv. 
    - The `perturbations` arg will be completely removed.
    - The meaning is "predict expression in `cell_type` at `time_point` if `perturbation` were set to `expression_level_after_perturbation`". However, methods may ignore some columns. For all aim 2 benchmarks, columns `cell_type`, `timepoint` will be placeholders and will be ignored. For GEARS, `expression_level_after_perturbation` is ignored. 

#### GGRN call to backends

- `predictions_metadata` will be passed through, except beforehand, `prediction_timescale` will be cartesian-producted with it, yielding one more column called `predictions_metadata`. To maintain DRY, any representation of `prediction_timescale` outside `predictions_metadata` should be ignored by individual backends or even better not passed to them. 