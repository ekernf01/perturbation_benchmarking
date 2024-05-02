## March-May 2024 updates to the benchmarking framework for timeseries prediction

- Want to be able to make predictions for any combo of `cell_type`, `timepoint`, `perturbation`, `expression_level_after_perturbation`, `prediction_timescale`. 
- For existing experiments, default is one prediction per test sample, with different `expression_level_after_perturbation`, `perturbation` but same `cell_type`, `timepoint`, `prediction_timescale`
- For TS, default will be only one `expression_level_after_perturbation` but all `cell_type`, `timepoint`, `perturbation`, `prediction_timescale`
- The interface between GGRN and benchmarking code will have to change -- cannot currently carry all this info. 

### Desired new interface: 

#### Pereggrn call to ggrn

- in GRN.predict(), the `perturbations` arg will be replaced by a new arg called `predictions_metadata`. It will be a pandas dataframe with columns `cell_type`, `timepoint`, `perturbation`, and `expression_level_after_perturbation`. 
    - It will default to `predictions.obs` or `starting_expression.obs` if those are provided.
    - It will be provided to Dockerized methods as a csv. 
    - The `perturbations` arg will be completely removed.
    - The meaning is "predict expression in `cell_type` at `time_point` if `perturbation` were set to `expression_level_after_perturbation`". However, methods may ignore some columns. For all aim 2 benchmarks, columns `cell_type`, `timepoint` will be placeholders and will be ignored. For GEARS, `expression_level_after_perturbation` is ignored. 

#### GGRN call to backends

- `predictions_metadata` will be passed through, except beforehand, `prediction_timescale` will be cartesian-producted with it, yielding one more column called `predictions_metadata`. To maintain DRY, any representation of `prediction_timescale` outside `predictions_metadata` should be ignored by individual backends (or even better, not passed to them). 

### Evaluation

#### Quantities

In an ideal scenario, each method would make specific and literal predictions of the form, "Starting from cell type $A$ and timepoint $T_1$ in the training data, we alter the dosage of gene $G$ to a value of $E$, and then we run the model forward a number of steps $S$, yielding a prediction expression vector of $X$ at time $T_2$." Realistically, this will not happen for a few reasons. 

- Modeling: 
    - Among the models we test, only PRESCIENT's timescale can be interpreted so straightforwardly. scKINETICS and CellOracle's timescales are not calibrated, meaning it is unclear how long a single forward step of the model takes.
    - Some models, notably CellOracle and PRESCIENT, are not meant to be interpreted gene by gene. They advertise only the ability to predict low-dimensional summaries of cell state, like cluster labels or PCA embeddings. 
- Data:
    - Timeseries and perturbation datasets usually have batch effects, making the quantitative expression levels not comparable. Sometimes the source of the effect is known: for example, in the definitive endoderm data, the timeseries data are from methanol-fixed cells and the perturb-seq from live cells. Other times, e.g. BETS, the source is not known but unsupervised analysis raises red flags.

Based on past work, here are some signals that we might have more success predicting.

- Short-term change in expression. Dictys uses this for model evaluation.
- Short-term change in PCA embeddings. CellOracle and scKINETICS focus on this.
- Short-term change in cell type proportions. CellOracle focuses on this.
- Long-term change in cell type proportions. PRESCIENT focuses on this.
- Cell type specific TF activity: rank the TF's by the magnitude of the change predicted upon knockout, e.g. cosine similarity of velocity before and after KO. scKINETICS and CellOracle focus on this.

When evaluating fold changes, the baseline expression has to be different between predicted and test data. For the test data, it should be a **test-set** control sample from the same timepoint and cell type. For the predictions, it should also be a **training-set** control sample from the same timepoint and cell type. 

#### Shapes and sizes

Given the new implementation and the lack of clarity about timescales, the shapes of the training data and test data are no longer guaranteed to match. For example, we may run CellOracle and make predictions at 1,2,3, and 10 time-steps, yielding 4 observations per starting state per perturbation. The test data may contain many replicates per perturbation. 

The default behavior should be to:

- Compare predictions to test data within each cell type and timepoint, averaging together all test samples.
- Evaluate predictions after different numbers of time-steps separately, even if they are returned inside the same AnnData object. 

Special behavior should be allowed for scenarios where there is a block in differentiation, as in the SMAD2 knockdown endoderm data. A typical successful evaluation might show the largest TF activity or backwards velocity not in the final cell state that is never reached, but rather in some earlier intermediate. I do not know how to automate this right now, and I might resort to manually evaluating summary plots of the predictions.