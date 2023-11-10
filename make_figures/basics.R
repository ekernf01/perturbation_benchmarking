library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")

unit_scale = function(x){
  x = x - min(x, na.rm = T)
  if(max(x, na.rm = T)>0){
    x = x / max(x, na.rm = T)
  }
  return(x)
}

collect_experiments = function(experiments, stratify_by_pert = T){
  X <- list()
  for (experiment in experiments) {
    if (stratify_by_pert){
      filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerPert.parquet")
    } else {
      filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerTarget.parquet")
    }
    X[[experiment]] <- arrow::read_parquet(filepath)
    X[[experiment]]$question %<>% as.character
    X[[experiment]]$refers_to %<>% as.character
    X[[experiment]]$color_by %<>% as.character
    X[[experiment]][["__index_level_0__"]] = NULL
    X[[experiment]][["__index_level_1__"]] = NULL
  }
  X <- bind_rows(X)
  return(X)
}


make_the_usual_labels_nice = function(X){
  try({X$factor_varied %<>% factor(levels = c("matching_method", "regression_method", "network_datasets"))}, silent = T)
  try({X$perturbation_dataset %<>% gsub("_", " ", .) %>% sub(" ", "\n", .) %>% gsub("γ", "g" , .)}, silent = T)
  try({X$network_datasets %<>% gsub("0$", "", .)})
  X$x = ifelse(X$factor_varied=="regression_method", X$regression_method, X$network_datasets)
  X$timescale_handling = paste0(X$matching_method, "_(", X$prediction_timescale, "-step)")
  X$x = ifelse(X$factor_varied=="matching_method", X$timescale_handling, X$x)
  X$x %<>% gsub("_", " ", .)
  the_usual_levels = gtools::mixedsort(unique(X$x))
  X %<>% mutate(x = factor(x, levels = unique(c("empty", "dense", "median", "mean", the_usual_levels))))
  return(X)
}

main_experiments = c("1.0_1",   "1.0_2",   "1.0_3",   "1.0_5",   "1.0_6",   "1.0_7",   "1.0_8",   "1.0_9",   "1.0_10",
                     "1.2.2_1", "1.2.2_2", "1.2.2_3", "1.2.2_5", "1.2.2_6", "1.2.2_7", "1.2.2_8", "1.2.2_9", "1.2.2_10",
                     "1.4.3_1", "1.4.3_2", "1.4.3_3", "1.4.3_5", "1.4.3_6", "1.4.3_7", "1.4.3_8", "1.4.3_9", "1.4.3_10")
# Panel A: initial benchmarks
{
  X = collect_experiments(main_experiments) 
  X %<>% make_the_usual_labels_nice
  X %<>% subset(x!="QuantileRegressor") # We never got quantile regression to run fast enough :(
  X %<>% subset(type_of_split=="interventional") 
  X %<>%
    group_by(x, perturbation_dataset, factor_varied) %>%
    summarise(across(starts_with("mae"), mean))
  ggplot(X) +
    geom_boxplot(aes(x = x, y = mae)) + 
    facet_grid(perturbation_dataset~factor_varied, scales = "free") + 
    labs(x = "", y = "Mean absolute error") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_y_continuous(labels = scales::label_number_si()) +
    geom_hline(data = subset(X, x %in% c("mean", "median", "empty", "dense")), aes(yintercept=mae, color = x)) 
  ggsave('plots/fig_basics_a.pdf', width = 6, height = 10)
}

# Same figure for cell type prediction
{
  X = collect_experiments(main_experiments) 
  X %<>% make_the_usual_labels_nice
  X %<>% subset(x!="QuantileRegressor") # We never got quantile regression to run fast enough :(
  X %<>% subset(type_of_split=="interventional") 
  X %<>%
    group_by(x, perturbation_dataset, factor_varied) %>%
    summarise(across(starts_with("cell"), mean))
  ggplot(X) +
    geom_boxplot(aes(x = x, y = cell_type_correct)) + 
    facet_grid(perturbation_dataset~factor_varied, scales = "free") + 
    labs(x = "", y = "Cell type accuracy") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_y_continuous(labels = scales::label_number_si()) + 
    geom_hline(data = subset(X, x %in% c("mean", "median", "empty", "dense")), aes(yintercept=cell_type_correct, color = x))
  ggsave('plots/fig_basics_a_ct.pdf', width = 6, height = 10)
}

# Panel b: simulations
{
  X = collect_experiments(c("1.9_0", "1.9_1", "1.9_2", "1.9_3", "1.9_4"))
  X %<>%
    group_by(network_datasets, perturbation_dataset, factor_varied) %>%
    summarise(across(starts_with("MAE"), mean))
  X %<>%    
    dplyr::mutate(perturbation_dataset = gsub("MARA_FANTOM4", "MARA FANTOM4", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("cellnet_human_Hg1332", "cellnet human Hg1332", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("cellnet_human_Hugene", "cellnet human Hugene", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("gtex_rna", "gtex rna", perturbation_dataset)) %>%
    dplyr::mutate(perturbation_dataset = gsub("celloracle_human", "celloracle human", perturbation_dataset)) %>%
    tidyr::separate(perturbation_dataset, into = c(".",  "true_network", "number_of_steps", "noise_sd"), sep = "_") %>%
    make_the_usual_labels_nice %>%
    dplyr::mutate(true_network = gsub(".*=", "", true_network)) %>%
    dplyr::mutate(noise_sd = gsub(".*=", "", noise_sd)) %>%
    dplyr::mutate(number_of_steps = gsub(".*=", "", number_of_steps)) %>%
    dplyr::mutate(number_of_steps = paste0(number_of_steps, " steps")) 
  X$is_true_network = X$true_network == X$x
  ggplot(X) +
    geom_point(aes(x = x, y = mae, color = is_true_network)) + 
    scale_color_manual(values = c("FALSE"="black", "TRUE"="RED")) +
    labs(x = "Network source", y = "Mean absolute error") +
    facet_wrap(~true_network, scales = "free_y", ncol = 1) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  ggsave('plots/fig_basics_simulation.pdf', width = 3.5, height = 6)
}

# Supplement: breakdown by target gene
check_if_beats_baselines = function(mae, x){
  if(any(x=="median") & any(x=="mean")){
    return((round(mae, 4) < round(mae[x=="median"], 4)) & (round(mae, 4) < round(mae[x=="mean"], 4) ) )
  } else {
    return(NA)
  }
}
check_mae_reduction = function(mae, x){
  mae_reduction = (pmin(mae[x=="median"], mae[x=="mean"]) - mae) / pmin(mae[x=="median"], mae[x=="mean"])
}
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

# Supplement: additional performance metrics 
{
  X = collect_experiments(main_experiments) %>% make_the_usual_labels_nice
  X %<>% subset(x!="QuantileRegressor")
  metrics = c("spearman", "mse_top_20", "mse_top_100", "mse_top_200",
               "mse", "mae", "proportion_correct_direction", "cell_type_correct")
  metrics_where_bigger_is_better = c("spearman", "proportion_correct_direction", "cell_type_correct")
  X %<>%
    group_by(x, perturbation_dataset, factor_varied) %>%
    summarise(across(metrics, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(metrics), names_to = "metric")
  # These metrics are on totally different scales
  X %<>% 
    subset(x != "QuantileRegressor") %>%
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

# Supplement: testing different data splits
{
  X = collect_experiments("1.8.4_0") 
  X %<>% make_the_usual_labels_nice
  metrics = c("spearman", "mse_top_20", "mse_top_100", "mse_top_200",
              "mse", "mae", "proportion_correct_direction")
  metrics_where_bigger_is_better = c("spearman", "proportion_correct_direction")
  X %<>%
    group_by(x, type_of_split, data_split_seed, regression_method) %>%
    summarise(across(metrics, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(metrics), names_to = "metric")
  X %<>% 
    group_by(metric, type_of_split, data_split_seed) %>%
    mutate(value = value*ifelse(metric %in% metrics_where_bigger_is_better, 1, -1)) %>%
    mutate(metric = paste(metric, ifelse(metric %in% metrics_where_bigger_is_better, "", "(inverted)"))) %>%
    mutate(scaled_value = unit_scale(value), is_best = scaled_value==1)
  ggplot(X) + 
    ggtitle("Data splitting strategies") +
    geom_tile(aes(y = metric,
                  fill = scaled_value,
                  x = regression_method
                  )) + 
    labs(fill = "Scaled metric\n(higher is better)") + 
    geom_point(data = subset(X, is_best), aes(x = regression_method, y = metric, color = is_best)) + 
    scale_color_manual(values = c("red")) +
    facet_wrap(paste0("seed=", data_split_seed)~type_of_split) 
  ggsave('plots/fig_data_splitting.pdf', width = 8, height = 8)
}

