import os
import pandas as pd
import pereggrn.experimenter as experimenter
import pereggrn_networks
import pereggrn_perturbations
pereggrn_networks.set_grn_location("../network_collection/networks")
pereggrn_perturbations.set_data_path("../perturbation_data/perturbations")
all_active_experiments = []
for experiment in os.listdir("experiments"):
    try:
        all_active_experiments.append(pd.DataFrame(
            {
                k:experimenter.validate_metadata(experiment, input_folder = "experiments")[k]
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
        

