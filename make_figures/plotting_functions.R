#' Load evaluation results from a list of benchmarking experiments.
#'
#' @param experiments
#' @param stratify_by_pert If true (default), collect results at per-perturbation resolution. Otherwise, collect results at per-target-gene resolution.
#'
collect_experiments = function(experiments, stratify_by_pert = T){
  X <- list()
  for (experiment in experiments) {
    if (stratify_by_pert){
      filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerPert.parquet")
    } else {
      filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerTarget.parquet")
    }
    try({
      X[[experiment]] <- arrow::read_parquet(filepath)
      X[[experiment]]$question %<>% as.character
      X[[experiment]]$refers_to %<>% as.character
      X[[experiment]]$color_by %<>% as.character
      X[[experiment]][["__index_level_0__"]] = NULL
      X[[experiment]][["__index_level_1__"]] = NULL
    })
    keep = c()
    for(i in seq_along(X[[experiment]])){
      if(any(!is.na(X[[experiment]][[i]]))){
        keep = c(keep, i)
      }
    }
    X[[experiment]] = X[[experiment]][keep]
  }
  X <- bind_rows(X)
  return(X)
}

#' Make plot labels nice and put elements in the desired order left to right.
#' 
make_the_usual_labels_nice = function(X){
  try({colnames(X) %<>% gsub("cell_type_correct", "cell_label_accuracy", .)}, silent = T)
  try({X$factor_varied %<>% gsub("_", " ", .)}, silent = T)
  try({X$factor_varied %<>% factor(levels = c("regression method", "network datasets", "matching method"))}, silent = T)
  try({X$perturbation_dataset %<>% gsub("_", " ", .) %>% sub(" ", "\n", .) %>% gsub("γ", "g" , .)}, silent = T)
  try({X$perturbation_dataset %<>% gsub("nakatake\nsimulated scrna", "nakatake\nscrna\nsimulated", .)}, silent = T) # Fits tighter
  try({X$network_datasets %<>% gsub("0$", "", .)})
  X$x = ifelse(X$factor_varied=="regression method", X[["regression_method"]], X[["network_datasets"]])
  X[["timescale handling"]] = paste0(X[["matching_method"]], " (", X[["prediction_timescale"]], "-step)")
  X$x = ifelse(X$factor_varied=="matching method", X[["timescale handling"]], X$x)
  X$x %<>% gsub("_", " ", .)
  the_usual_levels = gtools::mixedsort(unique(X$x))
  X %<>% mutate(x = factor(x, levels = unique(c("empty", "dense", "median", "mean", the_usual_levels))))
  return(X)
}

#' Shift and scale vector to have min 0, max 1
#' 
unit_scale = function(x){
  x = x - min(x, na.rm = T)
  if(max(x, na.rm = T)>0){
    x = x / max(x, na.rm = T)
  }
  return(x)
}

#' Check if any method beats the mean and median baselines. Lower is better.
#' 
check_if_beats_baselines = function(mae, x){
  if(any(x=="median") & any(x=="mean")){
    return((round(mae, 4) < round(mae[x=="median"], 4)) & (round(mae, 4) < round(mae[x=="mean"], 4) 
    ) )
  } else {
    return(NA)
  }
}

#' Shift and scale vector to have min 0, max 1
#' 
percent_change_from_best = function(x){
  m = max(x, na.rm = T)
  x = -abs(100*((m - x) / m))
  x = x %>% pmax(-10)
  return(x)
}

#' Plot all of our metrics.
#'
#' @param X Data as if from collect_experiments.
#' @param facet1 @param facet2 variables to facet the plot by.
#'
#' This function returns a ggplot object. This function produces a rigid plot format:
#' the X axis will be X$x, y will be the performance, color will indicate the performance metric. 
#' 
plot_all_metrics = function(
    X,
    facet1 = "perturbation_dataset",
    facet2 = "factor_varied",
    compare_across_rows = FALSE,
    metrics = c("spearman", "mse_top_20", "mse_top_100", "mse_top_200",
                "mse", "mae", "proportion_correct_direction", "cell_label_accuracy"),
    colorscale = NULL,
    metrics_where_bigger_is_better = c("spearman", "proportion_correct_direction", "cell_label_accuracy", "proportion correct direction", "cell type correct")
){
  X[["facet1"]] = X[[facet1]]
  X[["facet2"]] = X[[facet2]]
  X %<>%
    group_by(x, facet1, facet2) %>%
    summarise(across(metrics, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(metrics), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gtools::mixedsort(unique(X[["metric"]])))
  plot = ggplot(X) +
    geom_bar(aes(x = x, y = value, fill = metric), position = "dodge", stat = "identity") + 
    facet_grid(facet1~facet2, scales = "free") + 
    scale_fill_manual( values = setNames( colorscale, metrics %>% gsub("_", " ", .) ) ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  + 
    xlab("") + ylab("")
  return (plot)
}

#' Plot all of our metrics in a heatmap, shifted and scaled so that best is 1 and worst is 0.
#'
#' @param X Data as if from collect_experiments.
#' @param facet1 @param facet2 variables to facet the plot by.
#' @param compare_across_rows If TRUE, then experiments with different values of facet2 are compared directly. 
#' This is only appropriate if the data split is exactly the same, as in our main benchmark results figure.
#'
#' This function returns a ggplot object. This function produces a rigid plot format:
#' the X axis will be X$x, y will be the name of the evaluation metric, fill will be the value of the
#' evaluation metric.
#'
heatmap_all_metrics = function(
    X,
    facet1 = "perturbation_dataset",
    facet2 = "factor_varied",
    compare_across_rows = FALSE,
    metrics = c("spearman", "mse_top_20", "mse_top_100", "mse_top_200",
                "mse", "mae", "proportion_correct_direction", "cell_label_accuracy"),
    metrics_where_bigger_is_better = c("spearman", "proportion_correct_direction", "cell_label_accuracy", "proportion correct direction", "cell type correct")
){
  X[["facet1"]] = X[[facet1]]
  X[["facet2"]] = X[[facet2]]
  X %<>%
    group_by(x, facet1, facet2) %>%
    summarise(across(metrics, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(metrics), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  # Rescale metrics
  if (compare_across_rows) {
    X %<>% group_by(metric, facet1)
  } else {
    X %<>% group_by(metric, facet1, facet2)
  }
  X %<>%
    mutate(value = value*ifelse(metric %in% metrics_where_bigger_is_better, 1, -1)) %>%
    # mutate(metric = paste(metric, ifelse(metric %in% metrics_where_bigger_is_better, "", "(inverted)"))) %>%
    mutate(scaled_value = percent_change_from_best(value))
  X %<>% subset(!is.na(scaled_value))
  # plawt
  g = ggplot(X) +
    geom_tile(aes(x = x, y = metric, fill = scaled_value)) + 
    # scale_fill_gradient( breaks=c(0,1),labels=c("min","max"), limits=c(0,1)) +
    # geom_point(data = subset(X, is_best), aes(x = x, y = metric, color = is_best)) + 
    # scale_color_manual(values = c("red")) +
    # labs(x = "", y = "", fill = "Scaled\nvalue", color = "Best\nperformer") +
    labs(x = "", y = "", fill = "Percent change\nfrom best \n(capped at 10%)", color = "Best\nperformer") +
    
    facet_grid(facet1~facet2, scales = "free") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  return(g)
}

#' Plot all of our metrics in a heatmap, shifted and scaled so that best is 1 and worst is 0.
#'
#' @param X Data as if from collect_experiments.
#' @param facet1 @param facet2 variables to facet the plot by.
#' @param compare_across_rows If TRUE, then experiments with different values of facet2 are compared directly. 
#' This is only appropriate if the data split is exactly the same, as in our main benchmark results figure.
#'
#' This function returns a ggplot object. This function produces a rigid plot format:
#' the X axis will be X$x, y will be the name of the evaluation metric, fill will be the value of the
#' evaluation metric.
#'
plot_one_metric = function(
    X,
    facet1 = "perturbation_dataset",
    facet2 = "factor_varied",
    compare_across_rows = FALSE,
    metric = "mse", 
    threshold_outliers_at = Inf
){
  X[["facet1"]] = X[[facet1]]
  X[["facet2"]] = X[[facet2]]
  X %<>%
    group_by(x, facet1, facet2) %>%
    summarise(across(metric, mean))
  X = X[ X[[metric]] < threshold_outliers_at, ]
  g = ggplot(X) +
    geom_point(aes_string(x = "x", y = metric)) + 
    geom_hline(data = subset(X, x %in% c("mean", "median", "empty", "dense")), aes_string(yintercept=metric, color = "x")) +
    labs(color = "Non-informative\nbaselines") + 
    xlab("") +
    facet_grid(facet1~facet2, scales = "free") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  return(g)
}
