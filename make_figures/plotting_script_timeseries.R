library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
library(rjson)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")

# {
  X = collect_experiments(paste0("1.2.2_", 14:22)) %>% make_the_usual_labels_nice
  X %<>% subset(is_timescale_strict)
  heatmap_all_metrics(X,
                      metrics = c("spearman", "mse_top_20", "mse_top_100", "mse_top_200",
                                  "mse", "mae", "proportion_correct_direction", "cell_label_accuracy", 
                                  "distance_in_pca"),
                      compare_across_rows = T) 
  ggsave('timeseries_plots/matching.svg', width = 4, height = 8)
}