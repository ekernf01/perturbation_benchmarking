library("ggplot2")
library("dplyr")
library("stringr")
library("arrow")
library("magrittr")
library("rjson")
library("hdf5r")
# # Versions of AnnDataR and hdf5r we used:
# remotes::install_version("hdf5r", "1.3.11")
# devtools::install_github("scverse/anndataR", ref = "dbc6897")

setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")

# 1.2.2_14 is endoderm and it's complicated enough that the default evaluation code doesn't make biological sense. 
{
  X = collect_experiments(paste0("1.2.2_", c(14))) 
  X %<>% make_the_usual_labels_nice
  train = anndataR::read_h5ad(paste0("../../perturbation_data/perturbations/definitive_endoderm/train.h5ad"))
  train$obs[c("viz1", "viz2")] = train$obsm[["X_pca"]][,1:2]
  test  = anndataR::read_h5ad(paste0("../../perturbation_data/perturbations/definitive_endoderm/test.h5ad"))
  test$obs[c("viz1", "viz2")] = read.csv(paste0("../experiments/1.2.2_14/outputs/test_data_projection/0.csv"), row.names = 1)
  assertthat::are_equal(test$obs_names, rownames(test$obsm[["X_pca"]]))
  predictions  = anndataR::read_h5ad(paste0("../experiments/1.2.2_14/outputs/predictions/0.h5ad"))
  predictions$obs[c("viz1", "viz2")] = read.csv(paste0("../experiments/1.2.2_14/outputs/predictions_projection/0.csv"), row.names = 1)
  assertthat::are_equal(predictions$obs_names, rownames(predictions$obsm[["X_pca"]]))

  ggplot() +
    ggtitle("Train, test, and predicted data") +
    geom_point(aes(x=viz1, y=viz2), data = train$obs, color = "gray") +
    geom_point(aes(x=viz1, y=viz2), data = test$obs, color = "black") +
    geom_point(aes(x=viz1, y=viz2), data = predictions$obs, color = "red") + 
    theme_minimal()
  
  test$obs %<>% mutate(perturbation_simple = ifelse(perturbation %in% c("FOXH1", "SOX17", "SMAD2", "SMAD4"), as.character(perturbation), "other"))
  predictions$obs %<>% mutate(perturbation_simple = ifelse(perturbation %in% c("FOXH1", "SOX17", "SMAD2", "SMAD4"), as.character(perturbation), "other"))
  ggplot() +
    ggtitle("Test data distribution with train data backdrop") + 
    geom_point(aes(x=viz1, y=viz2), data = train$obs, color = "gray") +
    stat_bin_hex(aes(x=viz1, y=viz2, color = perturbation_simple, size = ..density..), geom = "point", data = test$obs) + 
    theme_minimal() 
  
  predictions$obs %>%
    merge(
      predictions$obs %>%
        subset(
          perturbation=="Scramble", 
          select = c("viz1", "viz2", "cell_type", "perturbation_type", "prediction_timescale"),
        ) %>% 
        rename(control_viz1 = viz1, control_viz2 = viz2) %>%
        unique,
      by = c("cell_type", "perturbation_type", "prediction_timescale"), 
      all.x = T, 
      all.y = F
    ) %>%
    ggplot() +
    ggtitle("One set of predictions (arrows point from simulated controls to simulated treatments)") + 
    geom_point(aes(x=viz1, y=viz2), data = train$obs, color = "gray") +
    geom_segment(aes(x = control_viz1, y = control_viz2,
                     xend = viz1, yend = viz2, 
                     color = perturbation_simple),
                 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
  theme_minimal() 
    
  
  ggrepel::geom_label_repel(
    aes( x = viz1_test, y = viz2_test, label = perturbation_x), 
    data = X %>%
        subset(
          cell_type != "endoderm", 
          select = c("cell_type", "viz1_test", "viz2_test", "perturbation_x")
        ) %>% 
        group_by(cell_type, perturbation_x) %>%
        summarize(viz1_test = mean(viz1_test), viz2_test = mean(viz2_test))
    ) +
    # facet_wrap(~perturbation_x) +
    ggtitle(paste0("All ", X$perturbation_dataset[500], " predictions"))
}

# Numbers 16-21 are other datasets that are simpler to deal with. 
{
  X = collect_experiments(paste0("1.2.2_", 16:21)) %>% make_the_usual_labels_nice
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, matching_method) %>%
    summarise(across(DEFAULT_METRICS, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(DEFAULT_METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = DEFAULT_METRICS)
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



