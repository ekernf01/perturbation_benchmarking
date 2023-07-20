library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")

# Figure <autoregressive>
collect_experiments = function(experiments){
  X <- list()
  for (experiment in experiments) {
    print(experiment)
    filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerPert.parquet")
    X[[experiment]] <- arrow::read_parquet(filepath)
    X[[experiment]]$color_by %<>% as.character
    X[[experiment]]$refers_to %<>% as.character
    X[[experiment]]$question %<>% as.character
  }
  X <- bind_rows(X)
  return(X)
}

    ""
    "low_dimensional_training"
X = collect_experiments(c("1.2.2_1","1.2.2_2","1.2.2_3","1.2.2_5","1.2.2_6","1.2.2_7","1.2.2_8","1.2.2_9","1.2.2_10"))
X$regression_method %<>% gsub("0$", "", .)
the_usual_levels = unique(X$regression_method)
X$regression_method %<>% factor(levels = unique(c("empty", "dense", "median", "mean", "celloracle human", the_usual_levels)))

X %<>%
  group_by(regression_method, low_dimensional_structure, perturbation_dataset, low_dimensional_training) %>%
  summarise(across(ends_with("benefit"), mean))
for(metric in c("mae_benefit")){
  ggplot(X) + 
    geom_point(aes_string(x = "regression_method", 
                          y = metric,
                          color='low_dimensional_training'), position = position_dodge(width=0.3)) + 
    labs(x='', 
         y = "MAE improvement over baseline") +
    facet_wrap(~perturbation_dataset) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}
ggsave(paste0('plots/fig_gears_', metric, '.pdf'), width = 10, height = 5)



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