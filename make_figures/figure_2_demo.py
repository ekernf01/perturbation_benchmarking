import numpy as np
import pandas as pd
import scanpy as sc
import pereggrn_perturbations
pereggrn_perturbations.set_data_path('../../perturbation_data/perturbations')
import os 
import altair as alt

print(os.listdir("../experiments/1.0_1/outputs"))
X = pd.read_parquet("../experiments/1.0_1/outputs/evaluationPerPert.parquet")
print(X.query("gene=='GATA3' & condition==0").T.to_csv())

# Metrics:
# pearsonCorr,0.8628334731048072
# spearmanCorr,0.6939872813635597
# logFCNorm2,116.28345489501952
# mae,0.24543215334415436
# mse,1915.5604908063688
# spearman,0.6103465858719491
# proportion_correct_direction,0.8333679186553227
# mse_top_20,207.51050154170025
# mse_top_100,529.9619580687759
# mse_top_200,752.8383862192422
# cell_type_correct,0.0