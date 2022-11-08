import os
import pandas as pd
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
import experimenter
importlib.reload(experimenter)
pd.concat(
    [
        pd.DataFrame({
            k:experimenter.validate_metadata(experiment, permissive = True)[0][k]
            for k in ["nickname", "unique_id", "refers_to", "readme"]
        }, index = [experiment])
        for experiment in os.listdir("experiments") if "1." in experiment
    ]
).sort_values("unique_id").to_csv("all_experiments.tsv", sep = "\t")
