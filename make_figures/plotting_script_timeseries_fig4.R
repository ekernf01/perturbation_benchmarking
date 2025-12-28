# environment setup
{
  library("viridis")
  library("gtools")
  library("ggplot2")
  library("dplyr")
  library("stringr")
  library("arrow")
  library("magrittr")
  library("rjson")
  library("ggrepel")
  setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
  source("plotting_functions.R")
  METRICS = "pearson"
}

# =========================== Figure 5: cell type-specific networks ===============================
{
  X = collect_experiments(paste0("1.4.1_", 0:4))
  X[X[["pearson"]]!=0, ]
  X %<>% add_network_cell_type_metadata
  X %<>% make_the_usual_labels_nice
  X %<>% subset(prediction_timescale==10)
  X$pearson %<>% tidyr::replace_na(0)
  X$pearson_top_20 %<>% tidyr::replace_na(0)
  X$pearson_top_100 %<>% tidyr::replace_na(0)
  X$pearson_top_200 %<>% tidyr::replace_na(0)
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, network_datasets, network_source, network_tissue, network_pretty, network_cell_type_matches ) %>%
    summarise(across(METRICS, mean, na.rm = T))
  X %<>% tidyr::pivot_longer(cols = all_of(METRICS), names_to = "metric")
  X %<>% subset(network_source != "humanbase")
  X$source = X$network_source %>% 
    gsub("_tissue", "", . ) %>% 
    gsub("cellnet_human.*", "cellnet_human", .) %>% 
    gsub("cellnet_mouse.*", "cellnet_mouse", .) %>% 
    gsub("0.5", "(filtered for \nposterior probability >50%)", ., fixed = T) %>% 
    gsub("_", " ", .)
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gtools::mixedsort(unique(X[["metric"]])))
  X[["network_cell_type"]] = c("Different", "Matched")[X[["network_cell_type_matches"]]+1]
  X[X[["network_tissue"]] %in% c("all"), "network_cell_type"] = "Not cell type-specific"
  X[X[["network_tissue"]] %in% c("dense"), "network_cell_type"] = "Fully connected network"
  X[X[["network_tissue"]] %in% c("empty"), "network_cell_type"] = "Empty network"
  X[["network_cell_type"]] %<>% factor(levels = c(
    "Different",
    "Matched"   ,
    "Not cell type-specific" ,
    "Fully connected network",
    "Empty network"
  ))
  X[["is_mismatch"]] = X[["network_cell_type"]]=="Different"
  X[["network_pretty"]] %<>% factor(levels = gtools::mixedsort(unique(X[["network_pretty"]])))

  for(mymetric in METRICS){
    plot = X %>% 
      subset(metric == mymetric) %>%
      ggplot() + 
      geom_point(
        aes(
          x = interaction(network_cell_type, source), 
          y = value, 
          color = source,
          alpha = network_cell_type,
          shape = network_cell_type
        ), 
        position = "dodge"
      ) +
      scale_alpha_manual(
        values = c(
          "Different" = 0.50,
          "Matched"   = 1,
          "Not cell type-specific" = 1,
          "Empty network" = 1,
          "Fully connected network" = 1
        )
      ) +
      ggrepel::geom_label_repel(
        aes(
          x = interaction(network_cell_type_matches, source),
          y = value,
          label = network_pretty
        ), data = subset(X, network_pretty %in% c() )
      ) +
      labs(fill = "Network cell type matches \nperturbation data cell type?", color = "Network source") +
      facet_wrap(~perturbation_dataset, scales = "free") +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  +
      xlab("") +
      ylab(mymetric) + 
      labs(shape = "Network cell type versus\nevaluation cell type", alpha = "Network cell type versus\nevaluation cell type")
    print(plot)
    ggsave(paste0('timeseries_plots/cell_type_specific_', mymetric, '.pdf'), plot, width = 6, height = 6)
  }
}



