library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
library(rjson)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")

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
