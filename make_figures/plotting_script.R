library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
library(rjson)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")

# Main experiments; all performance metrics 
main_experiments = c(  "1.0_1",   "1.0_2",   "1.0_3",   "1.0_5",   "1.0_6",   "1.0_7",   "1.0_8",   "1.0_9",   "1.0_10",   "1.0_11",   "1.0_12",   "1.0_13",
                       "1.2.2_1", "1.2.2_2", "1.2.2_3", "1.2.2_5", "1.2.2_6", "1.2.2_7", "1.2.2_8", "1.2.2_9", "1.2.2_10", "1.2.2_11", "1.2.2_12", "1.2.2_13",
                       "1.4.3_1", "1.4.3_2", "1.4.3_3", "1.4.3_5", "1.4.3_6", "1.4.3_7", "1.4.3_8", "1.4.3_9", "1.4.3_10", "1.4.3_11", "1.4.3_12", "1.4.3_13" )

{
  X = collect_experiments(main_experiments)
  X %<>% make_the_usual_labels_nice
  X %>% 
    subset( x=="mean") %>% 
    ggplot() + 
    # geom_jitter(aes(y = spearman, x = perturbation_dataset), alpha = 0.2) + 
    geom_boxplot(aes(y = spearman, x = perturbation_dataset), alpha = 1,
                 color = "black",
                 fill = NA,
                 outlier.shape = 'o') + 
    geom_hline(yintercept = 0, color = "red") + 
    ylab("Spearman correlation between\npredicted and observed fold change") + 
    ggtitle("Performance of 'mean' baseline") + 
    theme_minimal()    +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black")
    )  + 
    xlab("")
  ggsave('plots/fig_basics_correlation.pdf', width = 5, height = 4) 
  X %<>% distill_results(    
    facet1 = "perturbation_dataset", 
    facet2 = "factor_varied", 
    compare_across_rows = T
  )
  # stratified by dataset
  heatmap_all_metrics(
    X, 
    facet1 = "perturbation_dataset", 
    facet2 = "factor_varied"
  )  + 
    geom_vline(aes(xintercept = x), data = data.frame(x = 3.5, facet2 = factor("regression method", levels = levels(X$factor_varied)))) + 
    geom_vline(aes(xintercept = x), data = data.frame(x = 2.5, facet2 = factor("network datasets", levels = levels(X$factor_varied))))
  
  ggsave('plots/fig_basics_supp.pdf', width = 9, height = 14)
  # not stratified by dataset
  X %>% 
    # group_by(x, factor_varied, metric) %>%
    # summarize(scaled_value = median(scaled_value)) %>%
    mutate(metric = factor(metric, levels = gsub("_", " ", METRICS))) %>%
    ggplot() +
    geom_boxplot(aes(x = x, y = pmin(pmax(scaled_value, -100), 100))) + 
    labs(
      x = "", 
      y = "Percent change vs 'mean' baseline\nOriented so higher is better\nCapped at ±100%",
    ) +
    geom_hline(aes(yintercept = 0), colour = "red") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black")
    ) + 
    facet_grid(metric~factor_varied, scales = "free")
  ggsave('plots/fig_basics.pdf', width = 8, height = 14)
}

# Supplemental tables: stratifying performance by target gene or by perturbed gene
{
  evaluationPerTarget = collect_experiments(main_experiments, stratify_by_pert = FALSE)
  aggregated <- evaluationPerTarget %>% 
    subset(type_of_split == "interventional") %>%
    make_the_usual_labels_nice %>% 
    magrittr::extract(c("mae", "factor_varied", "x", "perturbation_dataset", "means", "variances", "variances_norm", 
              grep("degree", colnames(evaluationPerTarget), value = T),
              "n_exons", "pLI")) %>%
    tidyr::pivot_longer(
      cols = !(mae | factor_varied | x | perturbation_dataset),
      names_to = "property_of_gene", 
      values_to = "value"
    ) %>% 
    dplyr::select(perturbation_dataset, factor_varied, x, mae, property_of_gene, value) %>%
    mutate(value = as.numeric(value)) %>%
    na.omit() %>%
    group_by(property_of_gene, factor_varied, x, perturbation_dataset) %>%
    mutate(quintile = as.integer(cut(rank(value), breaks = 5))) %>%
    group_by(property_of_gene, factor_varied, x, perturbation_dataset, quintile) %>%
    summarise(value = median(value), mae = mean(mae)) %>%
    ungroup() %>%
    group_by(perturbation_dataset, property_of_gene, quintile)  %>%
    dplyr::mutate(
      beats_baselines = check_if_beats_baselines(mae, x), 
      mae_relative_reduction = percent_change_from_best_baseline(mae, x), 
      best_method = best_method(mae, x)
    ) %>%
    subset(!(x %in% c("mean", "median")))
  
  aggregated %>% 
    subset(
      grepl("in-degree.*Hugene", property_of_gene, ignore.case = T) |
        !grepl("degree", property_of_gene)
    ) %>%
    mutate(property_of_gene = gsub("in-degree_cellnet_human_Hugene", "degree", property_of_gene)) %>% 
    mutate(property_of_gene = gsub("variances_norm", "dispersion", property_of_gene)) %>%
    mutate(property_of_gene = gsub("means", "mean", property_of_gene)) %>%
    ggplot() +
    geom_boxplot(aes(x = property_of_gene, 
                     y = pmax(-25, mae_relative_reduction), 
                     group = interaction(quintile, property_of_gene), 
                     color = as.character(quintile)),
                 position = "dodge") + 
    facet_grid(perturbation_dataset~factor_varied, scales = "free") + 
    scale_color_viridis_d() +
    ylab("Percent improvement versus baseline MAE \n(higher is better\nCapped at -25)") + 
    ggtitle("Model performance on subsets of target genes") + 
    theme(axis.text.x = element_text(color = "black", angle = 90, vjust = 0.5, hjust = 1)) + 
    geom_hline(yintercept=0, color = "red") + 
    labs(color = "Quintile")
  ggsave("plots/fig_basics_stratify_targets.pdf", width = 10, height = 10)
  dim(aggregated) %>% write.csv('plots/fig_basics_stratify_targets_num_groups.csv')
}

{
  # Breakdown by perturbed gene
  evaluationPerPert = collect_experiments(main_experiments, stratify_by_pert = T)
  aggregated = evaluationPerPert %>% 
    subset(type_of_split == "interventional") %>%
    make_the_usual_labels_nice %>% 
    magrittr::extract(c("mae", "factor_varied", "x", "perturbation_dataset",
                        grep("degree", colnames(evaluationPerTarget), value = T),
                        "n_exons", "pLI")) %>%
    tidyr::pivot_longer(
      cols = !(mae | factor_varied | x | perturbation_dataset),
      names_to = "property_of_gene", 
      values_to = "value"
    ) %>% 
    dplyr::select(perturbation_dataset, factor_varied, x, mae, property_of_gene, value) %>%
    mutate(value = as.numeric(value)) %>%
    na.omit() %>%
    group_by(property_of_gene, factor_varied, x, perturbation_dataset) %>%
    mutate(quintile = as.integer(cut(rank(value), breaks = 5))) %>%
    ungroup %>%
    group_by(property_of_gene, factor_varied, x, perturbation_dataset, quintile) %>%
    summarise(value = median(value), mae = mean(mae)) %>% 
    dplyr::mutate(
      beats_baselines = check_if_beats_baselines(mae, x), 
      mae_relative_reduction = percent_change_from_best_baseline(mae, x), 
      best_method = best_method(mae, x)
    ) %>%
    subset(!(x %in% c("mean", "median")))
  
  dim(aggregated) %>% write.csv('plots/fig_basics_stratify_perts_num_groups.csv')
  long_data <- aggregated %>%
    dplyr::mutate(
      beats_baselines = check_if_beats_baselines(mae, x), 
      mae_relative_reduction = percent_change_from_best_baseline(mae, x), 
    ) %>% 
    subset(beats_baselines) %>% #nothing passes this filter LOLOLOL
    arrange(-mae_relative_reduction)                         
  long_data %>% write.csv('plots/fig_basics_stratify_perts.csv')
}

# Simulations
{
  X = collect_experiments(c("1.9_0", "1.9_1", "1.9_2", "1.9_3", "1.9_4")) 
  X %<>%    
    dplyr::mutate(perturbation_dataset = gsub("MARA_FANTOM4", "MARA FANTOM4", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("cellnet_human_Hg1332", "cellnet human Hg1332", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("cellnet_human_Hugene", "cellnet human Hugene", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("gtex_rna", "gtex rna", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("celloracle_human", "celloracle human", perturbation_dataset)) %>%
    tidyr::separate(perturbation_dataset, into = c(".",  "perturbation_dataset", "number_of_steps", "noise_sd"), sep = "_") %>%
    make_the_usual_labels_nice %>%
    dplyr::mutate(perturbation_dataset = gsub(".*=", "", perturbation_dataset)) %>%
    dplyr::mutate(noise_sd = gsub(".*=", "", noise_sd)) %>%
    dplyr::mutate(number_of_steps = gsub(".*=", "", number_of_steps)) %>%
    dplyr::mutate(number_of_steps = paste0(number_of_steps, " steps")) %>%
    subset(perturbation_dataset != "MARA\nFANTOM4") 
  X$is_true_network = X$perturbation_dataset == X$x
  X$num_observations_in_group = 1
  X$cell_label_accuracy = NA
  X$data_split_seed %<>% paste0("Data split: ", .)
  g = X %>%
    distill_results( facet2 = "data_split_seed",  compare_across_rows = F, baseline_condition = "dense" ) %>%
    heatmap_all_metrics(
      facet2 = "data_split_seed", 
      scales = "fixed", 
      cap = 5, 
      baseline_condition = "dense"
    )
  g +
    geom_vline(data = g$data %>% subset(x==gsub("\\s", " ", facet1)), color = "red", aes(xintercept = x) ) + 
    coord_fixed()
  ggsave('plots/fig_basics_simulation.pdf', width = 10, height = 8)
}

# Published methods
{
  X = collect_experiments(
    experiments = 
      c(
        "1.6.1_1",
        "1.6.1_2",
        "1.6.1_3",
        "1.6.1_4",
        "1.6.1_6",
        "1.6.1_7",
        "1.6.1_8",
        "1.6.1_9",
        "1.6.1_10",
        "1.6.1_11",
        "1.6.1_12",
        "1.6.1_13",
        "1.6.1_14",
        "1.6.1_15",
        "1.6.1_16",
        "1.6.1_17"
      )
  )
  X %<>% make_the_usual_labels_nice()
  X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
  method_tidy = c(
    "RidgeCV" = "GeneFormer", 
    "GEARS" = "GEARS",
    "DCDFG-spectral_radius-mlplr-False"="DCD-FG" ,
    "DCDFG-spectral_radius-linearlr-False"="NOTEARS-LR",
    "median"  = "median",                              
    "mean"  = "mean"   
  ) %>% rev
  X$regression_method = method_tidy[X$regression_method] %>% factor(levels = method_tidy)
  bulk_datasets = c(
    "nakatake", 
    "joung",
    "replogle2",
    "replogle2\ntf only",
    "replogle2\n large effect",
    "replogle3", 
    "replogle4", 
    "freimer",
    "frangieh\nIFNg v3"
  )
  no_raw_counts = c("dixit", "adamson", "norman")
  X %<>% subset(!((perturbation_dataset %in% bulk_datasets) & regression_method == "GEARS"))
  X %<>% subset(!((perturbation_dataset %in% no_raw_counts) & regression_method == "GeneFormer"))
  X$x = X$regression_method
  dir.create("plots", showWarnings = FALSE)
  X %>% 
    subset(starting_expression=="control") %>%
    distill_results(
      facet2 = "perturbation_dataset", 
      compare_across_rows = F 
    ) %>%
    heatmap_all_metrics(
      facet2 = "perturbation_dataset", 
      scales = "fixed", 
      do_wrap = T
    ) + coord_fixed()
  ggsave(filename = paste0("plots/fig_all_published.pdf"), width = 7, height = 8)
  # Not grouped by dataset
  X %>% 
    subset(starting_expression=="control") %>%
    distill_results(
      facet2 = "perturbation_dataset", 
      compare_across_rows = F 
    ) %>%
    # group_by(x, metric) %>%
    # summarize(scaled_value = median(scaled_value)) %>%
    mutate(metric = factor(metric, levels = gsub("_", " ", METRICS))) %>%
    ggplot() +
    geom_boxplot(aes(
      x = x,
      y = pmin(pmax(scaled_value, -100), 100))) + 
    labs(
      x = "", 
      y = "Percent change vs 'mean' baseline\nOriented so higher is better\nCapped at ±100%", 
    ) +
    geom_hline(aes(yintercept = 0), color = "red") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black")
    ) + 
    scale_fill_gradient2() +
    facet_wrap(~metric)
  ggsave(filename = paste0("plots/fig_all_published_simple.pdf"), width = 6, height = 5)
  
  X %>% 
    subset(regression_method!="GEARS") %>%
    subset(regression_method!="GeneFormer") %>%
    subset(unique_id!="1.6.1_17") %>%
    subset(unique_id!="1.6.1_11") %>%
    subset(unique_id!="1.6.1_8") %>%
    subset(unique_id!="1.6.1_9") %>%
    distill_results(
      facet1 = "starting_expression",
      facet2 = "perturbation_dataset", 
      compare_across_rows = F 
    ) %>%
    heatmap_all_metrics(
      ., 
      facet1 = "starting_expression",
      facet2 = "perturbation_dataset", 
      scales = "fixed", 
      do_wrap = F
    ) + 
    coord_fixed() 
  ggsave(filename = paste0("plots/fig_dcdfg.svg"), width = 12, height = 5)
}

# networks-only
{
  cell_type_matching = rjson::fromJSON(file = "../../accessory_data/cell_type_matching.json")
  do_subnetworks_match = function(perturbation_dataset, subnetwork){
    x = cell_type_matching$celltypes[perturbation_dataset]
    x %<>% unlist
    # If x is from a network that does not include any relevant subnetwork, the key might be missing from cell_type_matching$networks.
    # If x is from a network that does include a relevant subnetwork, this code runs as expected: check if this subnetwork is the relevant one.
    try(
      { 
        x = cell_type_matching$networks[x][[1]]
        # This incomprehensible one-liner converts a nested list of vectors to a tidy dataframe: {a:[1,2], b:2} becomes [[a,a,b], [1,2,2]].
        x = data.frame(network = Reduce(c, mapply(rep, names(x), sapply(x, length))), 
                       subnetwork = Reduce(c, x))
        x = paste(x$network, x$subnetwork)
        return(subnetwork %in% x)
      }, 
      silent = T
    )
    return(F)
  }
  X = collect_experiments(c("1.4.4_" %>% paste0(c(1:8)) ))
  X$cell_types_match = mapply(do_subnetworks_match, 
                              X$perturbation_dataset, 
                              X$network_datasets)
  X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
  X$perturbation_dataset %<>% gsub("γ", "g", .)
  X %<>% make_the_usual_labels_nice()
  my_levels = unique(c("frangieh\nIFNg v1", "frangieh\nIFNg v2", "frangieh\nIFNg v3", "nakatake", "nakatake\nscrna\nsimulated", X$perturbation_dataset))
  X$perturbation_dataset %<>% factor(levels = my_levels)
  X$network_datasets %<>% gsub(".parquet", "", .)
  X$network_datasets %<>% gsub(".csv_converted", "", .)
  X$network_datasets %<>% gsub("_top_filtered", "", .)
  X$network_source = X$network_datasets %>% 
    strsplit(" ") %>% 
    sapply(extract2, 1) 
  X$network_tissue = X$network_datasets %>% 
    paste("all") %>%
    strsplit(" ") %>%
    sapply(extract2, 2) %>% 
    tolower %>% 
    gsub("_", " ", .) %>%
    gsub("b lymphocyte", "bcell", .) %>%
    gsub(" memory", "", .) %>%
    gsub(" regulatory", "", .) %>%
    gsub(" conventional", "", .) %>%
    gsub(" naive", "", .) %>%
    gsub("retinal pigment epithelial", "rpe", .) %>%
    gsub("chronic lymphocytic leukemia", "", .) %>%
    gsub("chronic myelogenous leukemia", "", .) %>%
    gsub("suprapubic", "", .) %>%
    gsub("lower leg", "", .) %>%
    gsub("brain .*", "brain", .) %>%
    gsub("cell line", "", .) %>%
    gsub("muscleskel", "muscle", .) %>%
    gsub("pancreatic", "pancreas", .) 
  single_networks = c("celloracle_human",      "magnum_compendium_ppi" , "MARA_FANTOM4"     ,     "STRING"           ,     "magnum_compendium_32" )
  X$network_tissue[X$network_source %in% single_networks] = X$network_source[X$network_source %in% single_networks] 
  X$network_source[X$network_source %in% single_networks] = "other"
  X$network_pretty = paste(
    as.integer(as.factor(X$network_source)),
    X$network_tissue
  )
  for(dataset in X$perturbation_dataset %>% unique){
    current_X = subset(X, perturbation_dataset == dataset)
    networks_by_fc = current_X %>% 
      dplyr::group_by(network_pretty) %>%
      dplyr::summarise(fc_targets_vs_non_targets = median(fc_targets_vs_non_targets, na.rm = T)) %>%
      dplyr::arrange(fc_targets_vs_non_targets)
    current_X$network_pretty %<>% factor(levels = networks_by_fc$network_pretty)
    ggplot(current_X) + 
      geom_boxplot(outlier.shape = NA, 
                   aes(color = cell_types_match,
                       x = network_pretty, y = pmax(pmin(fc_targets_vs_non_targets, 0.5), -0.5))) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))    + 
      ggtitle("Perturbation response enrichment of known regulons") + 
      facet_wrap(~network_source, scales = "free", nrow = 1) + 
      geom_vline(xintercept = 0) + 
      ggtitle(dataset)
      ylab("Log fold change in target genes minus\nlog fold change in other genes") + 
      theme(axis.text.x = element_text(family = "mono", face = "bold"))
    ggsave(filename = paste0(paste0("plots/fig_network_only_", dataset, ".pdf")), width = 14, height = 8)
    networks_by_pvalue = current_X %>% 
      dplyr::group_by(network_pretty) %>%
      dplyr::summarise(fc_targets_vs_non_targets = median(-log10(pvalue_targets_vs_non_targets + 0.00001), na.rm = T)) %>%
      dplyr::arrange(fc_targets_vs_non_targets)
    current_X$network_pretty %<>% factor(levels = networks_by_pvalue$network_pretty)
    ggplot(current_X) + 
      geom_boxplot(outlier.shape = NA, 
                   aes(color = cell_types_match,
                       x = network_pretty, y = -log10(pvalue_targets_vs_non_targets))) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))    + 
      ggtitle("Perturbation response enrichment of known regulons") + 
      facet_wrap(~network_source, scales = "free", nrow = 1) + 
      geom_vline(xintercept = 0) + 
      ylab("-Log10 p-value of H0: \ntarget genes have same fc as non-targets") + 
      theme(axis.text.x = element_text(family = "mono", face = "bold")) + 
      ggtitle(dataset)
    ggsave(filename = paste0(paste0("plots/fig_network_only_", dataset, ".pdf")), width = 14, height = 8)
  }
}
