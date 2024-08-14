library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
library(rjson)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")

{
  X = collect_experiments(paste0("1.2.2_", 14:22)) %>% make_the_usual_labels_nice
  X$perturbation_dataset %<>% gsub("paul.", "paul", .)
  X %<>% subset(is_timescale_strict)
  plot_all_metrics(X,
                   metrics = c(
                     "overlap_top_20",
                     "overlap_top_100",
                     "overlap_top_200",
                     "pearson_top_20",
                     "pearson_top_100",
                     "pearson_top_200",
                     "proportion_correct_direction"
                   ),
                   colorscale = c(
                     "yellow", 
                     "goldenrod1",
                     "goldenrod4",
                     "cyan", 
                     "blue",
                     "navy", 
                     "black"
                   ),
                   compare_across_rows = T)
  ggsave('timeseries_plots/matching.svg', width = 4, height = 8)
}

