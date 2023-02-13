import os
import gc
import json
import yaml
import gc
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
from itertools import product
# Deal with various file paths specific to this project
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'network_collection', 'load_networks'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
import evaluator
import ggrn
import load_networks
import load_perturbations
importlib.reload(evaluator)
importlib.reload(load_networks)
importlib.reload(load_perturbations)
os.environ["GRN_PATH"]           = PROJECT_PATH + "network_collection/networks"
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbation_data/perturbations"

# Parse experiment metadata
def validate_metadata(
    experiment_name, 
    permissive = False
):
    with open(os.path.join("experiments", experiment_name, "metadata.json")) as f:
        metadata = json.load(f)
    if (not permissive) and ("is_active" in metadata.keys()) and (not metadata["is_active"]):
        raise ValueError("This experiment is marked as inactive. If you really want to run it, edit its metadata.json.")
    print("\n\nRaw metadata for experiment " + experiment_name + ":\n")
    print(yaml.dump(metadata))

    # If metadata refers to another experiment, go find missing metadata there.
    if "refers_to" in metadata.keys():
        with open(os.path.join("experiments", metadata["refers_to"], "metadata.json")) as f:
            other_metadata = json.load(f)
            try:
                assert other_metadata["is_active"], "Referring to an inactive experiment is not allowed."
            except KeyError:
                pass
        for key in other_metadata.keys():
            if key not in metadata.keys():
                metadata[key] = other_metadata[key]
    else:
        metadata["refers_to"] = None

    # Set defaults (None often defers to downstream code)
    defaults = {
        "pruning_parameter": None, 
        "pruning_strategy": "none",
        "network_prior": "ignore",
        "desired_heldout_fraction": 0.5,
        "type_of_split": "interventional",
        "regression_method": "RidgeCV",
        "time_strategy": "steady_state",
        "starting_expression": "control",
        "kwargs": None,
        "data_split_seed": 0,
        "baseline_condition": 0,
        "num_genes": 10000,
        "only_tfs_are_regulators": False,
        "merge_replicates": False,
        "network_datasets":{"dense":{}},
        "skip_bad_runs": True,
    }
    for k in defaults:
        if not k in metadata:
            metadata[k] = defaults[k]

    # network handling is complex; add some default behavior to reduce metadata boilerplate
    for netName in metadata["network_datasets"].keys():
        if not "subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["subnets"] = ["all"]
        if not "do_aggregate_subnets" in metadata["network_datasets"][netName].keys():
            metadata["network_datasets"][netName]["do_aggregate_subnets"] = False
    
    # Check all keys
    desired_keys = (
        # Experiment info
        "unique_id",
        "nickname",
        "readme",
        "question",
        "is_active",
        "factor_varied",    
        "color_by",
        "facet_by",
        # Data and preprocessing
        "network_datasets",
        "perturbation_dataset",
        "merge_replicates",
        "desired_heldout_fraction",
        "type_of_split",
        # Modeling decisions
        "pruning_parameter", 
        "pruning_strategy",
        "network_prior",
        "regression_method",
        "time_strategy",
    )
    missing = [k for k in desired_keys if k not in metadata.keys()]
    assert len(missing)==0, f"Metadata is missing some required keys: {' '.join(missing)}"
    
    # Check a few of the values
    assert experiment_name == metadata["unique_id"], "Experiment is labeled right"
    if not permissive:
        assert metadata["perturbation_dataset"] in set(load_perturbations.load_perturbation_metadata().query("is_ready=='yes'")["name"]), "perturbation data exist as named"
        for netName in metadata["network_datasets"].keys():
            assert netName in set(load_networks.load_grn_metadata()["name"]).union({"dense", "empty"}) or "random" in netName, "Networks exist as named"
            assert "subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"
            assert "do_aggregate_subnets" in metadata["network_datasets"][netName].keys(), "Optional metadata fields filled correctly"

    print("\nFully parsed metadata:\n")
    print(yaml.dump(metadata))

    return metadata

def lay_out_runs(
  networks: dict, 
  metadata: dict,
) -> pd.DataFrame:
    """Lay out the specific runs done in this experiment.

    Args:
    networks (dict): dict with string keys and LightNetwork values
    outputs (str): folder name to save results in
    metadata (dict): metadata for this Experiment, from metadata.json. See this repo's global README.

    Returns:
        pd.DataFrame: metadata on the different conditions in this experiment

    """
    metadata = metadata.copy() # We're gonna mangle it. :)
    # This will fuck up the cartesian product below.
    del metadata["kwargs"]
    # See experimenter.get_networks() to see how the metadata.json turns into this
    metadata["network_datasets"] = list(networks.keys())
    # This is just too bulky to want in the csv
    del metadata["readme"]
    # This can't just be cartesian-producted like everything else. We'll add it back after the product.
    baseline_condition = metadata["baseline_condition"]
    try:
        baseline_condition = baseline_condition.copy()
    except AttributeError:
        pass
    del metadata["baseline_condition"]

    # product splits strings if you don't wrap each in a list.
    for k in metadata.keys():
        if type(metadata[k]) != list:
            metadata[k] = [metadata[k]]
    # Combos 
    experiments =  pd.DataFrame(
        [row for row in product(*metadata.values())], 
        columns=metadata.keys()
    )
    # The dense network is represented as an empty network to save space.
    # Recommended usage is to set network_prior="ignore", otherwise the empty network will be taken literally. 
    for i in experiments.index:
        experiments.loc[i, "network_prior"] = \
        "ignore" if experiments.loc[i, "network_datasets"] == "dense" else experiments.loc[i, "network_prior"]

    experiments.index.name = "condition"
    experiments["baseline_condition"] = baseline_condition
    return experiments
  
def do_one_run(
    experiments: pd.DataFrame, 
    i: int, 
    train_data: anndata.AnnData, 
    test_data: anndata.AnnData, 
    networks: dict, 
    outputs: str,
    metadata: dict, 
    ) -> anndata.AnnData:
    """Do one run (fit a GRN model and make predictions) as part of this experiment.

    Args:
        experiments (pd.DataFrame): Output of lay_out_runs
        i (int): A value from the experiments.index
        Other args: see help(lay_out_runs)

    Returns:
        anndata.AnnData: Predicted expression
    """
    grn = ggrn.GRN(
        train=train_data, 
        network                 = networks[experiments.loc[i,'network_datasets']],
        only_tfs_are_regulators = experiments.loc[i,'only_tfs_are_regulators'],
    )
    grn.extract_tf_activity(method = "tf_rna")
    grn.fit(
        method                               = experiments.loc[i,"regression_method"], 
        cell_type_labels                     = None,
        cell_type_sharing_strategy           = "identical",
        network_prior                        = experiments.loc[i,"network_prior"],
        pruning_strategy                     = experiments.loc[i,"pruning_strategy"],
        pruning_parameter                    = experiments.loc[i,"pruning_parameter"],
        projection                           = "none",  
        time_strategy                        = experiments.loc[i,"time_strategy"],
        kwargs                               = metadata["kwargs"],
    )
    return grn

# TODO: move this to the network collection loader module?
def get_subnets(netName:str, subnets:list, target_genes = None, do_aggregate_subnets = False) -> dict:
    """Get gene regulatory networks for an experiment.

    Args:
        netName (str): Name of network to pull from collection, or "dense" or e.g. "random0.123" for random with density 12.3%. 
        subnets (list, optional): List of cell type- or tissue-specific subnetworks to include. 
        do_aggregate_subnets (bool, optional): If True, return has just one network named netName. If False,
            then returned dict has many separate networks named like netName + " " + subnet_name.

    Returns:
        dict: A dict containing base GRN's as LightNetwork objects (see the docs in the load_networks module in the networks collection.)
    """
    print("Getting network '" + netName + "'")
    gc.collect()
    if "random" in netName:
        networks = { 
            netName: load_networks.LightNetwork(
                df = evaluator.pivotNetworkWideToLong( 
                    load_networks.makeRandomNetwork( target_genes = target_genes, density = float( netName[6:] ) ) 
                ) 
            )
        }
    elif "empty" == netName or "dense" == netName:
        networks = { 
            netName: load_networks.LightNetwork(df=pd.DataFrame(index=[], columns=["regulator", "target", "weight"]))
        }
        if "dense"==netName:
            print("WARNING: for 'dense' network, returning an empty network. In GRN.fit(), use network_prior='ignore'. ")
    else:            
        networks = {}
        if do_aggregate_subnets:
            new_key = netName 
            if subnets[0]=="all":
                networks[new_key] = load_networks.LightNetwork(netName)
            else:
                networks[new_key] = load_networks.LightNetwork(netName, subnets)
        else:
            for subnet_name in subnets:
                new_key = netName + " " + subnet_name
                if subnets[0]=="all":
                    networks[new_key] = load_networks.LightNetwork(netName)
                else:
                    networks[new_key] = load_networks.LightNetwork(netName, [subnet_name])
    return networks

def filter_genes(expression_quantified: anndata.AnnData, num_genes: int, outputs: str) -> anndata.AnnData:
    """Filter a dataset, keeping only the top-ranked genes and the directly perturbed genes.
    The top N and perturbed genes may intersect, resulting in less than num_genes returned.
    For backwards compatibility with the DCD-FG benchmarks, we do not try to fix this.

    Args:
        expression_quantified (anndata.AnnData): _description_
        num_genes: Usually number. Expected non-numeric values are "all" or None or np.NaN, and for all those inputs, we keep all genes.

    Returns:
        anndata.AnnData: Input data, but maybe with fewer genes. 
    """
    assert "highly_variable_rank" in set(expression_quantified.var.columns)
    if num_genes is None or num_genes=="all" or np.isnan(num_genes):
        return expression_quantified

    # Perturbed genes
    targeted_genes = np.where(np.array(
        [1 if g in expression_quantified.uns["perturbed_and_measured_genes"] else 0
        for g in expression_quantified.var_names]
    ))[0]
    n_targeted = len(targeted_genes)

    # Top N minus # of perturbed genes
    try:      
        variable_genes = np.where(
            expression_quantified.var["highly_variable_rank"] < num_genes - n_targeted
        )[0]
    except:
        raise Exception("num_genes must act like a number w.r.t. < operator; received {num_genes}.")
    
    gene_indices = np.union1d(targeted_genes, variable_genes)
    gene_set = expression_quantified.var.index.values[gene_indices]
    pd.DataFrame({"genes_modeled": gene_set}).to_csv(os.path.join(outputs, "genes_modeled.csv"))
    return expression_quantified[:, gene_set].copy()


def set_up_data_networks_conditions(metadata, amount_to_do, outputs):
    """Set up the expression data, networks, and a sample sheet for this experiment."""
    # Data, networks, experiment sheet in that order because reasons
    perturbed_expression_data = load_perturbations.load_perturbation(metadata["perturbation_dataset"])
    try:
        perturbed_expression_data = perturbed_expression_data.to_memory()
    except ValueError: #Object is already in memory.
        pass
    try:
        perturbed_expression_data.X = perturbed_expression_data.X.toarray()
    except AttributeError: #Matrix is already dense.
        pass
    if metadata["merge_replicates"]:
        perturbed_expression_data = averageWithinPerturbation(ad=perturbed_expression_data)

    # Get networks
    networks = {}
    for netName in list(metadata["network_datasets"].keys()):
        networks = networks | get_subnets(
            netName, 
            subnets = metadata["network_datasets"][netName]["subnets"], 
            target_genes = perturbed_expression_data.var_names, 
            do_aggregate_subnets = metadata["network_datasets"][netName]["do_aggregate_subnets"]
        )

    # Lay out each set of params 
    experiments = lay_out_runs(
        networks=networks, 
        metadata=metadata,
    )
    try:
        old_experiments = pd.read_csv(os.path.join(outputs, "experiments.csv"), index_col=0)
        experiments.to_csv(        os.path.join(outputs, "new_experiments.csv") )
        experiments = pd.read_csv( os.path.join(outputs, "new_experiments.csv"), index_col=0 )
        if not experiments.equals(old_experiments):
            print(experiments)
            print(old_experiments)
            raise ValueError("Experiment layout has changed. Debug or delete previous experiments.csv. Saving new vs old for debugging.")
    except FileNotFoundError:
        pass
    experiments.to_csv( os.path.join(outputs, "experiments.csv") )

    # Simulate data if needed
    if "do_simulate" in metadata: 
        if amount_to_do=="evaluations":
            print("Finding previously simulated data.")
            perturbed_expression_data = sc.read_h5ad(os.path.join(outputs, "simulated_data.h5ad"))
        else:
            print("Simulating data.")
            grn = ggrn.GRN(
                train=perturbed_expression_data, 
                network=networks[metadata["do_simulate"]["network"]],
            )
            perturbed_expression_data = grn.simulate_data(
                [
                    (r[1][0], r[1][1]) 
                    for r in perturbed_expression_data.obs[["perturbation", "expression_level_after_perturbation"]].iterrows()
                ],
                effects = "uniform_on_provided_network",
                noise_sd = metadata["do_simulate"]["noise_sd"],
                seed = 0,
            )
            perturbed_expression_data.write_h5ad(os.path.join(outputs, "simulated_data.h5ad"))

    return perturbed_expression_data, networks, experiments

def splitDataWrapper(
    perturbed_expression_data,
    desired_heldout_fraction, 
    networks: dict, 
    network_behavior: str = "union", 
    type_of_split: str = "interventional" ,
    data_split_seed = None,
):
    """Split the data into train and test.

    Args:
        networks (dict): dict containing LightNetworks. Used to restrict what is allowed in the test set.
        network_behavior (str): How to restrict what is allowed in the test set.
    """
    if data_split_seed is None:
        data_split_seed = 0
    if network_behavior is None or network_behavior == "union":
        allowedRegulators = perturbed_expression_data.var_names
        if any([k not in {"dense", "empty"} for k in networks.keys()]):
            network_regulators = set.union(*[networks[key].get_all_regulators() for key in networks])
            allowedRegulators = allowedRegulators.intersection(network_regulators)
    else:
        raise ValueError(f"network_behavior currently only allows 'union'; got {network_behavior}")
    perturbed_expression_data_train, perturbed_expression_data_heldout = \
        splitData(
            perturbed_expression_data, 
            allowedRegulators, 
            desired_heldout_fraction = desired_heldout_fraction,
            type_of_split            = type_of_split,
            data_split_seed = data_split_seed,
        )
    return perturbed_expression_data_train, perturbed_expression_data_heldout

def splitData(adata, allowedRegulators, desired_heldout_fraction, type_of_split, data_split_seed):
    """Determine a train-test split satisfying constraints imposed by base networks and available data.
    
    A few factors complicate the training-test split. 

    - Perturbed genes may be absent from most base GRN's due to lack of motif information or ChIP data. 
        These perhaps should be excluded from the test data to avoid obvious failure cases.
    - Perturbed genes may not be measured. These perhaps should be excluded from test data because we can't
        reasonably separate their direct vs indirect effects.

    If type_of_split=="simple", we make no provision for dealing with the above concerns. The only restriction is that
    all controls go in the training set.
    If type_of_split=="interventional", the `allowedRegulators` arg can be specified in order to keep any user-specified
    problem cases out of the test data. No matter what, we still use those perturbed profiles as training data, hoping 
    they will provide useful info about attainable cell states and downstream causal effects. 

    For some collections of base networks, there are many factors ineligible for use as test data -- so many that 
    we use all the eligible ones for test and the only ineligible ones for training. 
    For other cases, such as dense base networks, we have more flexibility, so we send some perturbations to the 
    training set at random even if we would be able to use them in the test set.

    parameters:

    - adata (anndata.AnnData): Object satisfying the expectations outlined in the accompanying collection of perturbation data.
    - allowedRegulators (list or set): interventions allowed to be in the test set. 
    - type_of_split (str): if "interventional" (default), then any perturbation is placed in either the training or the test set, but not both. 
        If "simple", then we use a simple random split, and replicates of the same perturbation are allowed to go into different folds.

    """
    if data_split_seed is None:
        data_split_seed = 0
    # For a deterministic result when downsampling an iterable, setting a seed alone is not enough.
    # Must also avoid the use of sets. 
    if type_of_split == "interventional":
        get_unique_keep_order = lambda x: list(dict.fromkeys(x))
        allowedRegulators = [p for p in allowedRegulators if p in adata.uns["perturbed_and_measured_genes"]]
        testSetEligible   = [p for p in adata.obs["perturbation"] if     all(g in allowedRegulators for g in p.split(","))]
        testSetIneligible = [p for p in adata.obs["perturbation"] if not all(g in allowedRegulators for g in p.split(","))]
        allowedRegulators = get_unique_keep_order(allowedRegulators)
        testSetEligible   = get_unique_keep_order(testSetEligible)
        testSetIneligible = get_unique_keep_order(testSetIneligible)
        total_num_perts = len(testSetEligible) + len(testSetIneligible)
        eligible_heldout_fraction = len(testSetEligible)/(0.0+total_num_perts)
        if eligible_heldout_fraction < desired_heldout_fraction:
            print("Not enough profiles for the desired_heldout_fraction. Will use all available.")
            testSetPerturbations = testSetEligible
            trainingSetPerturbations = testSetIneligible
        elif eligible_heldout_fraction == desired_heldout_fraction: #nailed it
            testSetPerturbations = testSetEligible
            trainingSetPerturbations = testSetIneligible
        else:
            # Plenty of perts work for either.
            # Put some back in trainset to get the right size, even though we could use them in test set.
            numExcessTestEligible = int(np.ceil((eligible_heldout_fraction - desired_heldout_fraction)*total_num_perts))
            excessTestEligible = np.random.default_rng(seed=data_split_seed).choice(
                testSetEligible, 
                numExcessTestEligible, 
                replace = False)
            testSetPerturbations = [p for p in testSetEligible if p not in excessTestEligible]                      
            trainingSetPerturbations = list(testSetIneligible) + list(excessTestEligible) 
        # Now that the random part is done, we can start using sets. Order may change but content won't. 
        testSetPerturbations     = set(testSetPerturbations)
        trainingSetPerturbations = set(trainingSetPerturbations)
        adata_train    = adata[adata.obs["perturbation"].isin(trainingSetPerturbations),:]
        adata_heldout  = adata[adata.obs["perturbation"].isin(testSetPerturbations),    :]
        adata_train.uns[  "perturbed_and_measured_genes"]     = set(adata_train.uns[  "perturbed_and_measured_genes"]).intersection(trainingSetPerturbations)
        adata_heldout.uns["perturbed_and_measured_genes"]     = set(adata_heldout.uns["perturbed_and_measured_genes"]).intersection(testSetPerturbations)
        adata_train.uns[  "perturbed_but_not_measured_genes"] = set(adata_train.uns[  "perturbed_but_not_measured_genes"]).intersection(trainingSetPerturbations)
        adata_heldout.uns["perturbed_but_not_measured_genes"] = set(adata_heldout.uns["perturbed_but_not_measured_genes"]).intersection(testSetPerturbations)
        print("Test set size:")
        print(len(testSetPerturbations))
        print("Training set size:")
        print(len(trainingSetPerturbations))    
    elif type_of_split == "simple":
        np.random.seed(data_split_seed)
        train_obs = np.random.choice(
            replace=False, 
            a = adata.obs_names, 
            size = round(adata.shape[0]*(1-desired_heldout_fraction)), 
        )
        for o in adata.obs_names:
            if adata.obs.loc[o, "is_control"]:
                train_obs = np.append(train_obs, o)
        test_obs = [i for i in adata.obs_names if i not in train_obs]
        adata_train    = adata[train_obs,:]
        adata_heldout  = adata[test_obs,:]
        trainingSetPerturbations = set(  adata_train.obs["perturbation"].unique())
        testSetPerturbations     = set(adata_heldout.obs["perturbation"].unique())
        adata_train.uns[  "perturbed_and_measured_genes"]     = set(adata_train.uns[  "perturbed_and_measured_genes"]).intersection(trainingSetPerturbations)
        adata_heldout.uns["perturbed_and_measured_genes"]     = set(adata_heldout.uns["perturbed_and_measured_genes"]).intersection(testSetPerturbations)
        adata_train.uns[  "perturbed_but_not_measured_genes"] = set(adata_train.uns[  "perturbed_but_not_measured_genes"]).intersection(trainingSetPerturbations)
        adata_heldout.uns["perturbed_but_not_measured_genes"] = set(adata_heldout.uns["perturbed_but_not_measured_genes"]).intersection(testSetPerturbations)
    else:
        raise ValueError(f"`type_of_split` must be 'simple' or 'interventional'; got {type_of_split}.")
    return adata_train, adata_heldout


def averageWithinPerturbation(ad: anndata.AnnData, confounders = []):
    """Average the expression levels within each level of ad.obs["perturbation"].

    Args:
        ad (anndata.AnnData): Object conforming to the validity checks in the load_perturbations module.
    """
    if len(confounders) != 0:
        raise NotImplementedError("Haven't yet decided how to handle confounders when merging replicates.")

    perts = ad.obs["perturbation"].unique()
    new_ad = anndata.AnnData(
        X = np.zeros((len(perts), len(ad.var_names))),
        obs = pd.DataFrame(
            {"perturbation":perts}, 
            index = perts, 
            columns=ad.obs.columns.copy(),
        ),
        var = ad.var,
        dtype = np.float32
    )
    for p in perts:
        p_idx = ad.obs["perturbation"]==p
        new_ad[p,].X = ad[p_idx,:].X.mean(0)
        new_ad.obs.loc[p,:] = ad[p_idx,:].obs.iloc[0,:]
        try:
            new_ad.obs.loc[p,"expression_level_after_perturbation"] = ad.obs.loc[p_idx, "expression_level_after_perturbation"].mean()
        except:
            # If it's a multi-gene perturbation in the format "0,0,0", don't bother averaging
            # Hope to fix this eventually to average within each coord. 
            new_ad.obs.loc[p,"expression_level_after_perturbation"] = ad.obs.loc[p_idx, "expression_level_after_perturbation"][0]
    new_ad.obs = new_ad.obs.astype(dtype = {c:ad.obs.dtypes[c] for c in new_ad.obs.columns}, copy = True)
    new_ad.raw = ad.copy()
    new_ad.uns = ad.uns.copy()
    return new_ad


def downsample(adata: anndata.AnnData, proportion: float, seed = None, proportion_genes = 1):
    """Downsample training data to a given fraction, always keeping controls. 
    Args:
        adata (anndata.AnnData): _description_
        proportion (float): fraction of observations to keep. You may end up with a little extra because all controls are kept.
        proportion_genes (float): fraction of cells to keep. You may end up with a little extra because all perturbed genes are kept.
        seed (_type_, optional): RNG seed. Seed defaults to proportion so if you ask for 80% of cells, you get the same 80% every time.

    Returns:
        anndata.AnnData: Subsampled data.
    """
    if seed is None:
        seed = proportion
    np.random.seed(int(np.round(seed)))
    mask       = np.random.choice(a=[True, False], size=adata.obs.shape[0], p=[proportion,       1-proportion], replace = True)
    mask_genes = np.random.choice(a=[True, False], size=adata.var.shape[0], p=[proportion_genes, 1-proportion_genes], replace = True)
    adata = adata[adata.obs["is_control"] | mask, :].copy()
    perturbed_genes_remaining = set(adata.obs["perturbation"])
    adata = adata[:, [adata.var.index.isin(perturbed_genes_remaining)] | mask_genes].copy()
    print(adata.obs.shape)
    adata.uns["perturbed_but_not_measured_genes"] = set(adata.obs["perturbation"]).difference(  set(adata.var_names))
    adata.uns["perturbed_and_measured_genes"]     = set(adata.obs["perturbation"]).intersection(set(adata.var_names))
    return adata
