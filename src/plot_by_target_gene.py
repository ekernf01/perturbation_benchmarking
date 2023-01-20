import pandas as pd
import numpy as np
import os 
import sys
import altair as alt
PROJECT_PATH = '/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/'
import load_perturbations
os.environ["PERTURBATION_PATH"]  = PROJECT_PATH + "perturbation_data/perturbations"

sys.path.append(os.path.expanduser(os.path.join(PROJECT_PATH, 'perturbation_data', 'load_perturbations'))) 

evaluationPerTarget = pd.read_parquet("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/experiments/1.0_1/outputs/evaluationPerTarget.parquet")
