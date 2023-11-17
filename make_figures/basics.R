library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_code.R")
main_experiments = c("1.0_1",   "1.0_2",   "1.0_3",   "1.0_5",   "1.0_6",   "1.0_7",   "1.0_8",   "1.0_9",   "1.0_10",
                     "1.2.2_1", "1.2.2_2", "1.2.2_3", "1.2.2_5", "1.2.2_6", "1.2.2_7", "1.2.2_8", "1.2.2_9", "1.2.2_10",
                     "1.4.3_1", "1.4.3_2", "1.4.3_3", "1.4.3_5", "1.4.3_6", "1.4.3_7", "1.4.3_8", "1.4.3_9", "1.4.3_10")


# Panel a: all performance metrics 
{
  X = collect_experiments(main_experiments) %>% make_the_usual_labels_nice
  X$cell_type_correct = NA
  heatmap_all_metrics(X, compare_across_rows = T)
  ggsave('plots/fig_basics_metrics.pdf', width = 6, height = 8)
}

# Panel b: simulations
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
    dplyr::mutate(number_of_steps = paste0(number_of_steps, " steps")) 
  X$is_true_network = X$perturbation_dataset == X$x
  X$cell_type_correct = NA
  heatmap_all_metrics(X, compare_across_rows = T) 
  ggsave('plots/fig_basics_simulation.pdf', width = 3.5, height = 6)
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
        mae_relative_reduction = check_mae_reduction(mae, x), 
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
      mae_relative_reduction = check_mae_reduction(mae, x), 
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

# Supplement: testing different data splits
{
  X = collect_experiments("1.8.4_0") 
  X$x = X$regression_method
  heatmap_all_metrics(X, facet1 = "data_split_seed", facet2 = "type_of_split", compare_across_rows = FALSE)
  ggsave('plots/fig_data_splitting.pdf', width = 8, height = 8)
}

