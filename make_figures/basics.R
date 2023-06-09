library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")

collect_experiments = function(experiments){
  X <- list()
  for (experiment in experiments) {
    filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerPert.parquet")
    X[[experiment]] <- arrow::read_parquet(filepath)
    X[[experiment]]$question %<>% as.character
    X[[experiment]]$refers_to %<>% as.character
  }
  X <- bind_rows(X)
  return(X)
}
{
  X = collect_experiments(c("1.0_1", "1.0_2", "1.0_3", "1.0_5", "1.0_6", "1.4.3_1", "1.4.5_1", "1.4.7_1"))
  X$factor_varied %<>% factor(levels = c("regression_method", "network_datasets"))
  X$perturbation_dataset %<>% gsub("γ", "g" , .)
  X$x = ifelse(X$factor_varied=="regression_method", X$regression_method, X$network_datasets)
  X %<>%
    group_by(x, perturbation_dataset, factor_varied) %>%
    summarise(across(ends_with("benefit"), mean))
  
  ggplot(X) +
    geom_boxplot(aes(x = x, y = mae_benefit)) + 
    facet_grid(perturbation_dataset~factor_varied, scales = "free") + 
    labs(x = "", y = "Improvement in MAE over mean") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    geom_hline(aes(yintercept=0, color = "red")) + 
    theme(legend.position = "none")
  ggsave('plots/fig_basics_a.pdf', width = 5, height = 6)
}

{
  X = collect_experiments(c("1.9_0", "1.9_1", "1.9_2", "1.9_3"))
  X %<>%
    group_by(network_datasets, perturbation_dataset) %>%
    summarise(across(ends_with("benefit"), mean))
  X = X %>%
    tidyr::separate(perturbation_dataset, into = c(".",  "number_of_steps", "noise_sd"), sep = "_") %>%
    dplyr::mutate(number_of_steps = gsub(".*=", "", number_of_steps)) %>%
    dplyr::mutate(number_of_steps = paste0(number_of_steps, " steps")) %>%
    dplyr::mutate(noise_sd = gsub(".*=", "", noise_sd)) %>%
    dplyr::mutate(noise = ifelse(noise_sd=="1", "noise", "no noise")) 
  ggplot(X) +
    geom_boxplot(aes(x = network_datasets, y = mae_benefit)) + 
    labs(x = "Network source", y = "Improvement in MAE over empty network") +
    facet_wrap(number_of_steps~noise, scales = "free_y", ncol = 1) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    geom_hline(aes(yintercept=0, color = "red")) + 
    theme(legend.position = "none")
  ggsave('plots/fig_basics_simulation.pdf', width = 2.5, height = 6)
}
{
  evaluationPerTarget = arrow::read_parquet("../experiments/1.0_1/outputs/evaluationPerTarget.parquet")
  long_data <- evaluationPerTarget %>%
    tidyr::gather(var_name, value, c("means", "variances", "variances_norm")) %>%
    dplyr::select(regression_method, mae_benefit, property_of_gene = var_name, value) %>%
    mutate(value = as.numeric(value)) %>%
    na.omit() %>%
    group_by(property_of_gene, regression_method) %>%
    mutate(value_binned = cut(rank(value), breaks = 5)) %>%
    group_by(property_of_gene, regression_method, value_binned) %>%
    summarise(value = median(value), mae_benefit = median(mae_benefit)) %>%
    ungroup() %>%
    select(-value_binned)
  
  ggplot(long_data, aes(x = value, y = mae_benefit, color = regression_method)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ property_of_gene, ncol = 5, scales = "free_x") 
  ggsave('plots/fig_basics_stratify_targets.pdf', width = 8, height = 2.5)
}

{
  evaluationPerTarget = arrow::read_parquet("../experiments/1.0_1/outputs/evaluationPerPert.parquet")
  long_data <- evaluationPerTarget %>%
    tidyr::gather(var_name, value, c("logFC", "logFCNorm2", "pearsonCorr", "spearmanCorr")) %>%
    dplyr::select(regression_method, mae_benefit, property_of_gene = var_name, value) %>%
    mutate(value = as.numeric(value)) %>%
    na.omit() %>%
    group_by(property_of_gene, regression_method) %>%
    mutate(value_binned = cut(rank(value), breaks = 5)) %>%
    group_by(property_of_gene, regression_method, value_binned) %>%
    summarise(value = median(value), mae_benefit = median(mae_benefit)) %>%
    ungroup() %>%
    select(-value_binned)
  
  ggplot(long_data, aes(x = value, y = mae_benefit, color = regression_method)) +
    geom_point() +
    geom_line() +
    facet_wrap(~ property_of_gene, ncol = 5, scales = "free_x") 
  ggsave('plots/fig_basics_stratify_perts.pdf', width = 8, height = 2.5)
}
