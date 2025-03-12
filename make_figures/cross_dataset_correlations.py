import numpy as np
import pandas as pd
import scipy as sp
import itertools as it
import pereggrn_perturbations
import altair as alt

pereggrn_perturbations.set_data_path('../../perturbation_data/perturbations')

nakatake = pereggrn_perturbations.load_perturbation('nakatake')
joung = pereggrn_perturbations.load_perturbation('joung')
replogle2 = pereggrn_perturbations.load_perturbation('replogle2')
replogle3 = pereggrn_perturbations.load_perturbation('replogle3')
dixit = pereggrn_perturbations.load_perturbation('dixit')
adamson = pereggrn_perturbations.load_perturbation('adamson')

def cross_correlate_expression(data1, data2):

    # Subset data to genes common in both datasets
    common_genes = list(set(data1.var.index) & set(data2.var.index))
    data1 = data1[:, common_genes].copy()
    data2 = data2[:, common_genes].copy()

    # Compute baseline shared gene expressions
    ctrl1 = np.asarray(data1[data1.obs.is_control].X.mean(axis=0)).squeeze()
    ctrl2 = np.asarray(data2[data2.obs.is_control].X.mean(axis=0)).squeeze()
    ctrl1n = data1[data1.obs.is_control].obs.perturbation.unique()
    ctrl2n = data2[data2.obs.is_control].obs.perturbation.unique()
    
    # Focus on the shared genetic perturbations
    correlations = list()
    common_perts = list(set(data1.obs.perturbation) & set(data2.obs.perturbation) - set(ctrl1n) - set(ctrl2n))
    for p in common_perts:
        trt1 = np.asarray(data1[data1.obs.perturbation==p].X.mean(axis=0)).squeeze()
        trt2 = np.asarray(data2[data2.obs.perturbation==p].X.mean(axis=0)).squeeze()        
        lfc1 = trt1 - ctrl1 
        lfc2 = trt2 - ctrl2                # Log Fold Change - X is log-transformed
        correlations.append([
            sp.stats.pearsonr(lfc1, lfc2).statistic,
            sp.stats.spearmanr(lfc1, lfc2).statistic,
            p
        ])
    correlations = pd.DataFrame(correlations, columns=['Pearson', 'Spearman', 'Perturbation'])
    return correlations




CRISPRi = ['replogle2', 'replogle3', 'dixit', 'adamson']
CRISPRiCorrelations = list()
for d1, d2 in it.combinations(CRISPRi, r=2):
    corrs = cross_correlate_expression(eval(d1), eval(d2))
    corrs['Dataset 1'] = d1
    corrs['Dataset 2'] = d2
    CRISPRiCorrelations.append(corrs)
    print(d1, d2)
CRISPRiCorrelations = pd.concat(CRISPRiCorrelations)




OE = ['nakatake', 'joung']
OECorrelations = list()
for d1, d2 in it.combinations(OE, r=2):
    corrs = cross_correlate_expression(eval(d1), eval(d2))
    corrs['Dataset 1'] = d1
    corrs['Dataset 2'] = d2
    OECorrelations.append(corrs)
    print(d1, d2)
OECorrelations = pd.concat(OECorrelations)




CRISPRiCorrelationsLong = CRISPRiCorrelations.melt(id_vars=['Dataset 1', 'Dataset 2', 'Perturbation'],
                                                   value_vars=['Pearson', 'Spearman'],
                                                   var_name='Correlation Type',
                                                   value_name='Value')
OECorrelationsLong = OECorrelations.melt(id_vars=['Dataset 1', 'Dataset 2', 'Perturbation'],
                                         value_vars=['Pearson', 'Spearman'],
                                         var_name='Correlation Type',
                                         value_name='Value')
AllCorrelations = pd.concat([CRISPRiCorrelationsLong, OECorrelationsLong])
AllCorrelations["Datasets"] = AllCorrelations["Dataset 1"] + " vs " + AllCorrelations["Dataset 2"]
AllCorrelations = AllCorrelations.query("Datasets!='dixit vs adamson'")
alt.data_transformers.disable_max_rows()
chart = alt.Chart(
        AllCorrelations
    ).mark_boxplot(
    ).encode(
        x = "Datasets:N", 
        y = "Value:Q", 
        color = "Correlation Type:N",
        xOffset = "Correlation Type:N",
    ).properties(
        width=400, 
        height=200,
        title = "Cross-dataset correlations"
    ) 
chart.save("plots/cross_dataset_correlations.svg")
