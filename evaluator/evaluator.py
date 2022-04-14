"""evaluator.py is a collection of functions for making and testing predictions about expression fold change after genetic perturbations.
It dedicates particular attention to interfacing with CellOracle, a thought-provoking and flexible perturbation prediction method.
"""

def evaluateCausalModel(expression, predictedExpression, baseline, doPlots=False):
    plots = {}
    metrics = pd.DataFrame(index = predictedExpression.keys(), columns = ["spearman", "spearmanp"])
    for key in predictedExpression:
        observed = expression[expression.obs["perturbation"]==key,:].X.mean(axis=0) - baseline
        predicted = predictedExpression[key] - baseline
        if not any(np.isnan(predicted)):
            metrics.loc[key,["spearman","spearmanp"]] = spearmanr(observed, predicted)
        else:    
            metrics.loc[key,["spearman","spearmanp"]] = 0,1
        if doPlots:
            plots[key] = sns.scatterplot(x=observed, y=predicted)
            plots[key].set(title=key + " (Spearman rho="+ str(round(metrics.loc[key,"spearman"])) +")")
            plots[key].set_xlabel("Observed log fc", fontsize = 20)
            plots[key].set_ylabel("Predicted log fc", fontsize = 20)
            plots[key].plot()
            plt.figure()
    metrics["perturbation"] = metrics.index
    return metrics, plots            
    
    
def trainCausalModelAndPredict(expression, 
                               baseNetwork, 
                               memoizationName,
                               perturbations,
                               clusterColumnName,
                               pruningParameters = {"p":0.001, "threshold_number":2000}):
    """Train a causal model and predict outcomes of unseen perturbations."""
    print("Working on " + memoizationName)
    oracle = co.Oracle()
    oracle.import_anndata_as_raw_count(adata=ko_lab_esc_data,
                                   cluster_column_name=clusterColumnName,
                                   embedding_name="X_pca")
    oracle.import_TF_data(TF_info_matrix=baseNetwork)
    oracle.perform_PCA()
    n_comps = 50
    k = 1
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                          b_maxl=k*4, n_jobs=4)
    # Ridge regression pruning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=DeprecationWarning)
        links = oracle.get_links(cluster_name_for_GRN_unit=clusterColumnName, 
                                 alpha=10, 
                                 model_method = "bayesian_ridge",
                                 verbose_level=10,    
                                 test_mode=False, 
                                 n_jobs=14)
        links.filter_links(p=pruningParameters["p"], 
                           weight="coef_abs", 
                           threshold_number=pruningParameters["threshold_number"])
    
    links.links_dict = links.filtered_links.copy()
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
    oracle.to_hdf5(file_path=memoizationName)

    output = {}
    for goi, level in perturbations:
        print("Predicting " + goi)
        try:
            oracle.simulate_shift(perturb_condition={goi: level}, n_propagation=3, ignore_warning = True)
            output[goi] = oracle.adata["Control",:].layers['simulated_count'].squeeze()
        except ValueError as e:
            output[goi] = np.nan
            print("Prediction failed for " + goi + " with error " + str(e))
    return output
    
    
def networkEdgesToMatrix(networkEdges, regulatorColumn=0, targetColumn=1):
    """Reformat a network from a two-column dataframe to the way that celloracle needs its input."""
    X = pd.crosstab(networkEdges.iloc[:,targetColumn], networkEdges.iloc[:,regulatorColumn])
    X = 1.0*(X > 0)
    X = X.rename_axis('gene_short_name').reset_index()
    X = X.rename_axis('peak_id').reset_index()
    return X

def countMatrixEdges(network):
    """Count the number of connections in a network that is formatted how CO expects it to be formatted."""
    return 1.0*network.iloc[:,2:].sum().sum()

def makeRandomNetwork(density = 0, seed = 0):
    """Generate a random network formatted the way that celloracle needs its input."""
    np.random.seed(seed)
    return pd.DataFrame(
            np.random.binomial(
                n = 1, 
                p=density,
                size=(
                    len(targetGenes), 
                    len(humanTFs['HGNC symbol'])
                )
            ),
            columns = humanTFs['HGNC symbol'], 
            index = targetGenes
        ).rename_axis('gene_short_name').reset_index().rename_axis('peak_id').reset_index()

