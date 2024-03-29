library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
library(rjson)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")


# Main experiments; all performance metrics 
main_experiments = c("1.0_1",   "1.0_2",   "1.0_3",   "1.0_5",   "1.0_6",   "1.0_7",   "1.0_8",   "1.0_9",   "1.0_10",
                     "1.2.2_1", "1.2.2_2", "1.2.2_3", "1.2.2_5", "1.2.2_6", "1.2.2_7", "1.2.2_8", "1.2.2_9", "1.2.2_10",
                     "1.4.3_1", "1.4.3_2", "1.4.3_3", "1.4.3_5", "1.4.3_6", "1.4.3_7", "1.4.3_8", "1.4.3_9", "1.4.3_10")
{
  X = collect_experiments(main_experiments) %>% make_the_usual_labels_nice
  plot_one_metric(X, compare_across_rows = T) 
  heatmap_all_metrics(X, compare_across_rows = T) 
  ggsave('plots/fig_basics.pdf', width = 8, height = 10)
}

# Supplemental tables: stratifying performance by target gene or by perturbed gene
{
  evaluationPerTarget = collect_experiments(main_experiments, stratify_by_pert = FALSE)
  long_data <- evaluationPerTarget %>% 
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
    group_by(perturbation_dataset, property_of_gene, quintile) %>%
    dplyr::mutate(
      beats_baselines = check_if_beats_baselines(mae, x), 
        mae_relative_reduction = percent_change_from_best(mae, x), 
    ) %>% 
    subset(beats_baselines) %>% 
    arrange(-mae_relative_reduction)
  long_data %>% write.csv("plots/fig_basics_stratify_targets.csv")
  
  # Breakdown by perturbed gene
  evaluationPerPert = collect_experiments(main_experiments, stratify_by_pert = T)
  long_data <- evaluationPerPert %>% 
    subset(type_of_split == "interventional") %>%
    make_the_usual_labels_nice %>% 
    magrittr::extract(c("mae", "factor_varied", "x", "perturbation_dataset",
                        "logFC", "logFCNorm2", "pearsonCorr", "spearmanCorr",
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
    group_by(perturbation_dataset, property_of_gene, quintile) %>%
    dplyr::mutate(
      beats_baselines = check_if_beats_baselines(mae, x), 
      mae_relative_reduction = percent_change_from_best(mae, x), 
    ) %>% 
    subset(mae_relative_reduction > 0.05) %>% 
    arrange(-mae_relative_reduction)                         
  long_data %>% write.csv('plots/fig_basics_stratify_perts.csv')
  long_data %>% 
    extract(c("property_of_gene", "quintile")) %>% 
    table %>% 
    as.data.frame %>% 
    arrange(-Freq) %>% 
    subset(Freq>0) %>%
    write.csv('plots/fig_basics_stratify_perts_summary.csv')
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
  X$cell_type_correct = NA
  g = heatmap_all_metrics(X, facet2 = "data_split_seed", compare_across_rows = F) 
  g + geom_vline(data = g$data %>% subset(x==gsub("\\s", " ", facet1)), color = "green", aes(xintercept = x) )
  ggsave('plots/fig_basics_simulation.pdf', width = 8, height = 8)
}

# Runtime analysis
{
  conditions = read.csv("../experiments/5_0/outputs/conditions.csv")
  X = list.files("../experiments/5_0/outputs/train_resources", full.names = T) %>% 
    lapply(read.csv) %>% 
    data.table::rbindlist() %>%
    dplyr::rename(condition = X) %>%
    merge(conditions, by = "condition") %>%
    dplyr::mutate(method = ifelse(feature_extraction=="geneformer", feature_extraction, regression_method)) %>%
    dplyr::mutate(peak.RAM = gsub("KB","E3", peak.RAM)) %>%
    dplyr::mutate(peak.RAM = gsub("MB","E6", peak.RAM)) %>%
    dplyr::mutate(peak.RAM = gsub("GB","E9", peak.RAM)) %>%
    dplyr::mutate(peak.RAM = gsub("B","", peak.RAM))  %>%
    dplyr::mutate(peak.RAM = as.numeric(peak.RAM)) 
    
  ggplot(X) + 
    geom_point(aes(x = method, y = walltime..seconds., color = num_genes)) + 
    scale_y_log10() + 
    ylab("Walltime (seconds)") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle("Compute time on Nakatake with different numbers of genes predicted") 
  ggplot(X) + 
    geom_point(aes(x = method, y = peak.RAM, color = num_genes)) + 
    scale_y_log10() + 
    ylab("Peak RAM (bytes)") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
    ggtitle("Peak RAM consumption on Nakatake with different numbers of genes predicted") 
  ggsave("plots/fig_ram.pdf", width = 5, height = 5)
}

# Geneformer
{
  X = collect_experiments(c("1.3.3_1","1.3.3_2",
                            "1.3.3_3","1.3.3_5", "1.3.3_6", "1.3.3_7", "1.3.3_8", "1.3.3_9", "1.3.3_10"))
  X$regression_method %<>% gsub("0$", "", .)
  X$regression_method %<>% gsub("RidgeCV", "Regress on\nGeneFormer embeddings", .)
  X %<>% make_the_usual_labels_nice()
  heatmap_all_metrics(X, facet2 = "perturbation_dataset", facet1 = "factor_varied")
  ggsave(paste0('plots/fig_geneformer.pdf'), width = 10, height = 4)
}

# GEARS
{
  X = collect_experiments(c(
    "1.4.2_1",
    "1.4.2_2",
    "1.4.2_3",
    "1.4.2_4",
    "1.4.2_12",
    "1.4.2_13",
    "1.4.2_14"
  )
  )
  X %<>% make_the_usual_labels_nice()
  X$regression_method %<>% gsub("0$", "", .)
  the_usual_levels = unique(X$regression_method)
  X$regression_method %<>% factor(levels = unique(c("empty", "dense", "median", "mean", "celloracle human", the_usual_levels)))
  X$x = X$regression_method
  X$facet2 = with(X, paste0(desired_heldout_fraction*100, "% heldout\nseed=", data_split_seed))
  heatmap_all_metrics(X, facet2 = "perturbation_dataset", facet1 = "facet2")
  ggsave(paste0('plots/fig_gears.pdf'), width = 8, height = 4)
}

# DCD-FG
{
  X = collect_experiments(c("1.6.1_1", "1.6.1_3", "1.6.1_6",  "1.6.1_7", "1.6.1_16", "1.6.1_2",
                            "1.6.1_10", "1.6.1_11", "1.6.1_12", "1.6.1_13", "1.6.1_14", "1.6.1_15", "1.6.1_16"))
  X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
  method_tidy = c(
    "DCDFG-spectral_radius-mlplr-False"="DCD-FG" ,
    "median"  = "median",                              
    "DCDFG-spectral_radius-linearlr-False"="NOTEARS-LR",
    "mean"  = "mean"   
  ) 
  X$regression_method = method_tidy[X$regression_method] %>% factor(levels = c("median", "mean", "NOTEARS-LR", "DCD-FG"))
  X$perturbation_dataset %<>% gsub("γ", "g", .)
  X %<>% make_the_usual_labels_nice()
  my_levels = unique(c("frangieh\nIFNg v1", "frangieh\nIFNg v2", "frangieh\nIFNg v3", "nakatake", "nakatake\nscrna\nsimulated", X$perturbation_dataset))
  X$perturbation_dataset %<>% factor(levels = my_levels)
  X$x = X$regression_method
  heatmap_all_metrics(X, facet2 = "perturbation_dataset", facet1 = "starting_expression", compare_across_rows = F)
  dir.create("plots", showWarnings = FALSE)
  ggsave(filename = paste0("plots/fig_dcdfg.pdf"), width = 12, height = 4)
}

# Exact repeat of a DCD-FG experiment
{
  X = collect_experiments(c("1.6.1_15", "1.6.1_19"))
  X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
  method_tidy = c(
    "DCDFG-spectral_radius-mlplr-False"="DCD-FG" ,
    "median"  = "median",                              
    "DCDFG-spectral_radius-linearlr-False"="NOTEARS-LR",
    "mean"  = "mean"   
  ) 
  X$regression_method = method_tidy[X$regression_method] %>% factor(levels = c("median", "mean", "NOTEARS-LR", "DCD-FG"))
  X$perturbation_dataset %<>% gsub("γ", "g", .)
  X %<>% make_the_usual_labels_nice()
  my_levels = unique(c("frangieh\nIFNg v1", "frangieh\nIFNg v2", "frangieh\nIFNg v3", "nakatake", "nakatake\nscrna\nsimulated", X$perturbation_dataset))
  X$perturbation_dataset %<>% factor(levels = my_levels)
  X$x = X$regression_method
  heatmap_all_metrics(X, facet2 = "unique_id", facet1 = "starting_expression", compare_across_rows = F)
  ggsave(filename = paste0("plots/fig_dcdfg_repeat.pdf"), width = 6, height = 3)
}

{
  X = collect_experiments(c("1.6.1_17"))
  X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
  method_tidy = c(
    "DCDFG-spectral_radius-mlplr-False"="DCD-FG" ,
    "median"  = "median",                              
    "DCDFG-spectral_radius-linearlr-False"="NOTEARS-LR",
    "mean"  = "mean"   
  ) 
  X$regression_method = method_tidy[X$regression_method] %>% factor(levels = c("median", "mean", "NOTEARS-LR", "DCD-FG"))
  X$perturbation_dataset %<>% gsub("γ", "g", .)
  X %<>% make_the_usual_labels_nice()
  my_levels = unique(c("frangieh\nIFNg v1", "frangieh\nIFNg v2", "frangieh\nIFNg v3", "nakatake", "nakatake\nscrna\nsimulated", X$perturbation_dataset))
  X$perturbation_dataset %<>% factor(levels = my_levels)
  X$x = X$regression_method
  heatmap_all_metrics(X, facet1 = "perturbation_dataset", facet2 = "data_split_seed", compare_across_rows = F)
  dir.create("plots", showWarnings = FALSE)
  ggsave(filename = paste0("plots/fig_dcdfg_followup.pdf"), width = 6, height = 3)
}

{
  X = collect_experiments(c("1.6.1_18"))
  X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
  method_tidy = c(
    "DCDFG-spectral_radius-mlplr-False"="DCD-FG" ,
    "median"  = "median",                              
    "DCDFG-spectral_radius-linearlr-False"="NOTEARS-LR",
    "mean"  = "mean"   
  ) 
  X$regression_method = method_tidy[X$regression_method]
  X$regression_method %<>% paste(formatC(X$pruning_parameter, format = "e", digits = 0)) 
  X$regression_method %<>% factor(levels = rev(gtools::mixedsort(unique(X$regression_method))))
  X$perturbation_dataset %<>% gsub("γ", "g", .)
  X %<>% make_the_usual_labels_nice()
  my_levels = unique(c("frangieh\nIFNg v1", "frangieh\nIFNg v2", "frangieh\nIFNg v3", "nakatake", "nakatake\nscrna\nsimulated", X$perturbation_dataset))
  X$perturbation_dataset %<>% factor(levels = my_levels)
  X$x = X$regression_method
  heatmap_all_metrics(X, facet1 = "perturbation_dataset", facet2 = "data_split_seed", compare_across_rows = F)
  dir.create("plots", showWarnings = FALSE)
  ggsave(filename = paste0("plots/fig_dcdfg_tuning.pdf"), width = 6, height = 3)
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
