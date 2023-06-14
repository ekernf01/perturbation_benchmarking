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
    X[[experiment]]$color_by %<>% as.character
    
  }
  X <- bind_rows(X)
  return(X)
}

make_the_usual_labels_nice = function(X){
  try({X$factor_varied %<>% factor(levels = c("regression_method", "network_datasets"))}, silent = T)
  try({X$perturbation_dataset %<>% gsub("_", " ", .) %>% sub(" ", "\n", .) %>% gsub("γ", "g" , .)}, silent = T)
  try({X$network_datasets %<>% gsub("0$", "", .)})
  X$x = ifelse(X$factor_varied=="regression_method", X$regression_method, X$network_datasets)
  X$x %<>% gsub("_", " ", .)
  the_usual_levels = unique(X$x)
  X %<>% mutate(x = factor(x, levels = unique(c("empty", "dense", "median", "mean", "celloracle human", the_usual_levels))))
  return(X)
}
main_experiments = c("1.0_1",     "1.0_2",   "1.0_3",   "1.0_5",   "1.0_6",   "1.0_7",   "1.0_8",   "1.0_9", "1.0_10",
                     "1.4.3_1", "1.4.3_2", "1.4.3_3", "1.4.3_5", "1.4.3_6", "1.4.3_7", "1.4.3_8", "1.4.3_9", "1.4.3_10")
{
  X = collect_experiments(main_experiments) 
  X %<>% make_the_usual_labels_nice
  X %<>% subset(x!="QuantileRegressor")
  X %<>%
    group_by(x, perturbation_dataset, factor_varied) %>%
    summarise(across(ends_with("benefit"), mean))
  X$factor_varied %<>% factor(levels = c("regression_method", "network_datasets"))
  ggplot(X) +
    geom_boxplot(aes(x = x, y = mae_benefit)) + 
    facet_grid(perturbation_dataset~factor_varied, scales = "free") + 
    labs(x = "", y = "Improvement in MAE over mean") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    geom_hline(aes(yintercept=0, color = "red")) + 
    theme(legend.position = "none")
  ggsave('plots/fig_basics_a.pdf', width = 5, height = 8)
}

{
  X = collect_experiments(c("1.9_0", "1.9_1", "1.9_2", "1.9_3"))
  X %<>%
    group_by(network_datasets, perturbation_dataset, factor_varied) %>%
    summarise(across(ends_with("benefit"), mean))
  X = X %>%
    tidyr::separate(perturbation_dataset, into = c(".",  "number_of_steps", "noise_sd"), sep = "_") %>%
    make_the_usual_labels_nice %>%
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
    make_the_usual_labels_nice %>%
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
  ggsave('plots/fig_basics_stratify_targets.pdf', width = 10, height = 2)
}

{
  evaluationPerTarget = arrow::read_parquet("../experiments/1.0_1/outputs/evaluationPerPert.parquet")
  long_data <- evaluationPerTarget %>% 
    make_the_usual_labels_nice %>%
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
  ggsave('plots/fig_basics_stratify_perts.pdf', width = 10, height = 2)
}

{
  X = collect_experiments(main_experiments) %>% make_the_usual_labels_nice
  X %<>% subset(x!="QuantileRegressor")
  metrics = c("spearman", "mse_top_20", "mse_top_100", "mse_top_200",
               "mse", "mae", "proportion_correct_direction")
  metrics_where_bigger_is_better = c("spearman", "proportion_correct_direction")
  X %<>%
    group_by(x, perturbation_dataset, factor_varied) %>%
    summarise(across(metrics, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(metrics), names_to = "metric")
  unit_scale = function(x){
    x = x - min(x)
    if(max(x)>0){
      x = x / max(x)
    }
    return(x)
  }
  # These metrics are on totally different scales
  X %<>% 
    subset(x != "regression_metric") %>%
    group_by(metric, perturbation_dataset) %>%
    mutate(value = value*ifelse(metric %in% metrics_where_bigger_is_better, 1, -1)) %>%
    mutate(metric = paste(metric, ifelse(metric %in% metrics_where_bigger_is_better, "", "(inverted)"))) %>%
    mutate(scaled_value = unit_scale(value), is_best = scaled_value==1)
  ggplot(X) +
    geom_tile(aes(x = x, y = metric, fill = scaled_value)) + 
    geom_point(data = subset(X, is_best), aes(x = x, y = metric, color = is_best)) + 
    scale_color_manual(values = c("red")) +
    facet_grid(perturbation_dataset~factor_varied, scales = "free") + 
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  ggsave('plots/fig_basics_metrics.pdf', width = 6, height = 8)
}
