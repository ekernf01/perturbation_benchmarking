"""evaluator.py is a collection of functions for making and testing predictions about expression fold change after genetic perturbations.
It dedicates particular attention to interfacing with CellOracle, a thought-provoking and flexible perturbation prediction method.
"""
import celloracle as co
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import gc
from scipy.stats import spearmanr as spearmanr

def evaluateCausalModel(heldout, predictions, baseline, experiments, experiment_name, factor_varied, classifier = None):
    """Compile plots and tables comparing heldout data and predictions for same. Saves output to f"results/{experiment_name}".

    Args:
        heldout (AnnData): Test data in the format specified by the complementary collection of perturbation data.
        predictions (Dict): Keys are perturbation names; values are np arrays of predicted expression or np.nan if no prediction made. 
        baseline (np.Array): Control expression for use in calculating fold change. 
        classifier (sklearn.LogisticRegression): Optional, to judge results on cell type accuracy. 
        experiments (pd.DataFrame): Metadata for the different combinations used in this experiment. 
        factor_varied (String): Plots are automatically stratified based on this column of "experiments". 
        experiment_name (String): Saves output to f"results/{experiment_name}".
    """
    # Get spearman and classifier accuracy 
    evaluationResults = {}
    for i in predictions:
        evaluationResults[i] = \
            evaluateOnePrediction(heldout, 
                                  predictions,   
                                  baseline = baseline,     
                                  doPlots=False, 
                                  classifier=classifier)[0]
        evaluationResults[i]["index"] = i
    evaluationResults = pd.concat(evaluationResults)
    evaluationResults = evaluationResults.merge(experiments, how = "left")
    evaluationResults = pd.DataFrame(evaluationResults.to_dict())
    # Mark anything where predictions were unavailable
    noPredictionMade = evaluationResults.iloc[[x==0 for x in evaluationResults["spearman"]],:]['perturbation']
    noPredictionMade = set(noPredictionMade)
    noPredictionMade
    evaluationResults["somePredictionRefused"] = evaluationResults["perturbation"].isin(noPredictionMade) 
    # Save main results
    evaluationResults.to_csv("results/"+ experiment_name +"/evaluationResults.csv")
    stripchartMainFig = sns.stripplot(y = "spearman", 
                                                x=factor_varied, 
                                                data = evaluationResults)
    stripchartMainFig.set(ylabel="Spearman correlation \npredicted log fold change vs observed")
    stripchartMainFig.figure.savefig(f'results/{EXPERIMENT_NAME}/stripchart.pdf')
    plt.show()
    for readout in "spearman", "cell_fate_correct":
        meanSEPlot = sns.pointplot(factor_varied, y=readout, data=evaluationResults, dodge=True, join=False)
        meanSEPlot.set(title="Mean and SE of " + readout)
        meanSEPlot.figure.savefig(f'results/{EXPERIMENT_NAME}/MeanSEPlot{readout}.pdf')
        plt.show()
        plt.figure()
    # Which genes are best/worst?
    hardest = evaluationResults[evaluationResults["spearman"].idxmax(),"perturbation"]
    easiest = evaluationResults[evaluationResults["spearman"].idxmin(),"perturbation"]
    evaluationResults.loc[evaluationResults["perturbation"]==hardest,:].to_csv("results/"+ experiment_name +"/hardest.csv")
    evaluationResults.loc[evaluationResults["perturbation"]==easiest,:].to_csv("results/"+ experiment_name +"/easiest.csv")
    return

def evaluateOnePrediction(expression, predictedExpression, baseline, doPlots=False, classifier = None):
    '''Compare observed against predicted, for expression, fold-change, or cell type.

            Parameters:
                    expression (dict of 1d numpy arrays): 
                        the actual expression (log-scale). 
                    predictedExpression (dict of 1d numpy arrays): 
                        the cellOracle prediction (log-scale). Elements of the dict may be np.nan for 
                        missing predictions.
                    baseline (numpy array): 
                        control expression level (log-scale)
                    classifier (sklearn logistic regression classifier): 
                        optional machine learning classifier to assign cell fate. 
                        Must have a predict() method capable of taking a value from expression or predictedExpression and returning a single class label. 

            Returns:
                    Pandas DataFrame with Spearman correlation between predicted and observed 
                    log fold change over control.
    '''
    "log fold change using Spearman correlation and (optionally) cell fate classification."""
    plots = {}
    metrics = pd.DataFrame(index = predictedExpression.keys(), columns = ["spearman", "spearmanp", "cell_fate_correct"])
    for key in predictedExpression:
        observed = expression[expression.obs["perturbation"]==key,:].X.mean(axis=0)
        predicted = predictedExpression[key]
        if type(predicted) is float and np.isnan(predicted):
            metrics.loc[key,["spearman","spearmanp", "cell_fate_correct"]] = 0,1,np.nan
        else:
            metrics.loc[key,["spearman","spearmanp"]] = spearmanr(observed - baseline, predicted - baseline)
            if classifier is not None:
                class_observed = classifier.predict(np.reshape(observed, (1, -1)))[0]
                class_predicted = classifier.predict(np.reshape(predicted, (1, -1)))[0]
                metrics.loc[key,"cell_fate_correct"] = 1.0*(class_observed==class_predicted)            
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
                               perturbations,
                               clusterColumnName,
                               memoizationName=None,
                               pruningParameters = {"p":0.001, "threshold_number":2000}):
    """Train a causal model and predict outcomes of unseen perturbations.

    Args:
        expression (AnnData): AnnData object; training data as described in this project's collection of perturbation data.
        baseNetwork (pd.DataFrame): Base GRN in the format expected by CellOracle.
        perturbations (list): List of tuples [(gene, expression_level)] where the gene is fixed at the expression level.
        clusterColumnName (_type_): Categorical column from AnnData input to use as cluster labels. 
        memoizationName (_type_, optional): _description_. Defaults to None.
        pruningParameters (dict, optional): _description_. Defaults to {"p":0.001, "threshold_number":2000}.

    Returns:
        _type_: _description_
    """
    sc.pl.umap(expression, color = clusterColumnName)
    if memoizationName:
        print("Working on " + memoizationName)
    # Default output if something goes wrong
    output = {}
    for goi, level in perturbations:
        output[goi] = np.nan
        
    # Memoization
    try:         
        oracle = co.load_hdf5(file_path=memoizationName)
        print("Memoized results found.")
    except (ValueError, AttributeError) as e:
        print("Memoized results not found with error " + str(e))
        
        # Object setup
        oracle = co.Oracle()
        oracle.import_anndata_as_raw_count(adata=expression,
                                       cluster_column_name=clusterColumnName,
                                       embedding_name="X_pca")
        baseNetwork = makeNetworkDense(baseNetwork)        
        oracle.import_TF_data(TF_info_matrix=baseNetwork)
        oracle.perform_PCA()
        n_comps = 50
        k = 1
        oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                              b_maxl=k*4, n_jobs=4)
        
        # Training
        try:
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
            if memoizationName:
                oracle.to_hdf5(file_path=memoizationName)
        except Exception as e: #LinAlgError: SVD did not converge?
            print("Training failed with error " + str(e))
            return output
    
    # Prediction
    print("Predicting " + goi)
    for goi, level in perturbations:
        print(goi, end=' ')
        try:
            oracle.simulate_shift(perturb_condition={goi: level}, n_propagation=3, ignore_warning = True)
            output[goi] = oracle.adata["Control",:].layers['simulated_count'].squeeze()
        except ValueError as e:
            print("\nPrediction failed for " + goi + " with error " + str(e))
    
    # Free memory & return
    del oracle
    gc.collect()
    return output
    
    
def networkEdgesToMatrix(networkEdges, regulatorColumn=0, targetColumn=1):
    """Reformat a network from a two-column dataframe to the way that celloracle needs its input."""
    X = pd.crosstab(networkEdges.iloc[:,targetColumn], networkEdges.iloc[:,regulatorColumn])
    del networkEdges
    gc.collect()
    X = 1.0*(X > 0)
    X = X.rename_axis('gene_short_name').reset_index()
    X = X.rename_axis('peak_id').reset_index()
    X = makeNetworkSparse(X, 0.0)
    gc.collect()
    return X

def countMatrixEdges(network):
    """Count the number of connections in a network that is formatted how CO expects it to be formatted."""
    return 1.0*network.iloc[:,2:].sum().sum()

humanTFs = pd.read_csv("../accessory_data/humanTFs.csv")
targetGenes = co.data.load_human_promoter_base_GRN()["gene_short_name"]

def makeRandomNetwork(density = 0, seed = 0, TFs = humanTFs['HGNC symbol'], targetGenes = targetGenes ):
    """Generate a random network formatted the way that celloracle needs its input."""
    np.random.seed(seed)
    X = pd.DataFrame(
            np.random.binomial(
                n = 1, 
                p=density,
                size=(
                    len(targetGenes), 
                    len(TFs)
                )
            ),
            columns = TFs, 
            index = targetGenes
        )
    X.rename_axis('gene_short_name', inplace=True)
    X.reset_index(inplace=True)
    X.rename_axis('peak_id', inplace=True)
    X.reset_index(inplace=True)
    # CellOracle's preferred format wastes gobs of memory unless you sparsify.
    X = makeNetworkSparse(X, round(density))
    gc.collect()
    return X

def makeNetworkSparse(X, defaultValue):
    """Save memory by making a sparse representation of a base network"""
    X.iloc[:,2:] = X.iloc[:,2:].astype(pd.SparseDtype("float", defaultValue))
    return X

def makeNetworkDense(X):
    """Undo makeNetworkSparse"""
    X.iloc[:, 2:] = np.array(X.iloc[:, 2:])   #undo sparse representation         
    return X

def splitData(adata, allowedRegulators, minTestSetSize = 250, perturbationColName = "perturbation"):
    """Determine a train-test split satisfying constraints imposed by base networks and available data.
    
A few factors complicate the training-test split. 

- Perturbed genes may be absent from most base GRN's due to lack of motif information or ChIP data. These are excluded from the test data to avoid obvious failure cases.
- Perturbed genes may not be measured. These are excluded from test data because we don't know to what extent they were overexpressed.

In both cases, we still use those perturbed profiles as training data, hoping they will provide useful info about attainable cell states and downstream causal effects. 

For some collections of base networks, there are many factors ineligible for use as test data -- so many that we use all the eligible ones for test and the only ineligible ones for training. For other cases, such as dense base networks, we have more flexibility, so we send some perturbations to the training set at random even if we would be able to use them in the test set.

parameters:
    - adata: AnnData object satisfying the expectations outlined in the accompanying collection of perturbation data.
    - allowedRegulators: list or set of features allowed to be in the test set. In CellOracle, this is usually constrained by motif/chip availability. 

    """
    allowedRegulators = set(allowedRegulators)
    allowedRegulators = allowedRegulators.difference(adata.uns["perturbed_but_not_measured_genes"])
    if len(allowedRegulators) < minTestSetSize:
        raise ValueError(f"minTestSetSize too big; set to {minTestSetSize} but {len(allowedRegulators)} available.")
    testSetPerturbations     = set(adata.obs[perturbationColName]).intersection(allowedRegulators)
    trainingSetPerturbations = set(adata.obs[perturbationColName]).difference(allowedRegulators)
    if len(trainingSetPerturbations) < minTestSetSize:
        swap = np.random.default_rng(seed=0).choice(list(testSetPerturbations), 
                                               minTestSetSize - len(trainingSetPerturbations), 
                                               replace = False)
        testSetPerturbations = testSetPerturbations.difference(swap)
        trainingSetPerturbations = trainingSetPerturbations.union(swap)
    adata_heldout  = adata[adata.obs[perturbationColName].isin(testSetPerturbations),    :]
    adata_train    = adata[adata.obs[perturbationColName].isin(trainingSetPerturbations),:]
    adata_train.obs['perturbation'].unique()
    perturbationsToPredict = [(gene, adata_heldout[sample, gene].X[0,0]) for sample,gene in 
                              enumerate(adata_heldout.obs[perturbationColName])] 
    print("Example perturbations formatted as \n (gene, expression after perturbation)")
    print(perturbationsToPredict[0:5])
    print("Test set size:")
    print(len(testSetPerturbations))
    print("Training set size:")
    print(len(trainingSetPerturbations))
    return adata_train, adata_heldout, perturbationsToPredict
