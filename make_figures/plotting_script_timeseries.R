library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
library(rjson)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")
DEFAULT_METRICS = c(
  "overlap_top_20",
  "overlap_top_100",
  "overlap_top_200",
  "pearson_top_20",
  "pearson_top_100",
  "pearson_top_200",
  "proportion_correct_direction"
)
{
  X = collect_experiments(paste0("1.2.2_", 14:22)) %>% make_the_usual_labels_nice
  X$perturbation_dataset %<>% gsub("paul.", "paul", .)
  X %<>% subset(is_timescale_strict)
  
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, matching_method) %>%
    summarise(across(DEFAULT_METRICS, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(DEFAULT_METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gtools::mixedsort(unique(X[["metric"]])))
  X[["prediction_timescale"]] %<>% as.character()
  X[["prediction_timescale"]] %<>% factor(levels = gtools::mixedsort(unique(X[["prediction_timescale"]])))
  plot = ggplot(X) +
    geom_bar(aes(x = prediction_timescale, y = value, fill = matching_method), position = "dodge", stat = "identity") +
    facet_grid(metric~perturbation_dataset, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
    ylab("")
  print(plot)
  ggsave('timeseries_plots/matching.svg', plot, width = 4, height = 8)
}

{
  X = collect_experiments(paste0("1.4.1_", 0:4))
  X %<>% subset(is_timescale_strict)
  X[X[["pearson_top_20"]]!=0, ]
  X %<>% add_network_cell_type_metadata
  X %<>% make_the_usual_labels_nice
  X %<>% subset(is_timescale_strict)
  X %<>% subset(prediction_timescale==10)
  X$pearson %<>% tidyr::replace_na(0)
  X$pearson_top_20 %<>% tidyr::replace_na(0)
  X$pearson_top_100 %<>% tidyr::replace_na(0)
  X$pearson_top_200 %<>% tidyr::replace_na(0)
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, network_datasets, network_source, network_tissue, network_pretty, network_cell_type_matches ) %>%
    summarise(across(DEFAULT_METRICS, mean, na.rm = T))
  X %<>% tidyr::pivot_longer(cols = all_of(DEFAULT_METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gtools::mixedsort(unique(X[["metric"]])))
  X %<>% subset(network_source != "humanbase")
  X$source = X$network_source %>% 
    gsub("_tissue", "", . ) %>% 
    gsub("cellnet_human.*", "cellnet_human", .) %>% 
    gsub("cellnet_mouse.*", "cellnet_mouse", .)
  for(mymetric in c("overlap top 20", "pearson top 20")){
    plot = X %>% 
      subset(metric == mymetric) %>%
      ggplot() + 
      geom_bar(
        aes(
          x = network_pretty, 
          y = value, 
          color = source, 
          fill = c("Different", "Matched")[network_cell_type_matches+1]
        ), 
        stat = "identity", 
        position = "dodge") +
      scale_fill_manual(values = c("gray", "black")) +
      labs(fill = "Network cell type matches \nperturbation data cell type?", color = "Network source") +
      facet_wrap(~perturbation_dataset, scales = "free", ncol = 1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
      xlab("") + 
      ylab(mymetric)
    print(plot)
    ggsave(paste0('timeseries_plots/cell_type_specific_', mymetric, '.svg'), plot, width = 12, height = 12)
  }
}



