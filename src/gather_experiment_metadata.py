import os
import pandas as pd
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
os.chdir(PROJECT_PATH + "perturbation_benchmarking")
import importlib
import sys
sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_benchmarking', 'src'))) 
import experimenter
importlib.reload(experimenter)
all_active_experiments = []
for experiment in os.listdir("experiments"):
    try:
        all_active_experiments.append(pd.DataFrame(
            {
                k:experimenter.validate_metadata(experiment, permissive = True)[k]
                for k in ["nickname", "refers_to", "readme"]
            }, 
            index = [experiment]
        ))
    except:
        all_active_experiments.append(pd.DataFrame(
            {
                k:"Could not validate the metadata -- likely an inactive experiment."
                for k in ["nickname", "refers_to", "readme"]
            }, 
            index = [experiment]
        ))
pd.concat(all_active_experiments).sort_index().to_csv("all_experiments.tsv", sep = "\t", index = True)
print("Done. See results in all_experiments.tsv.")
        

