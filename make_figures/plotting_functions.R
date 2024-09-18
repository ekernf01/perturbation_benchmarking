# Ordered by perturbation type, then by duration. See table 1. 
DATASET_ORDER = c(
  "nakatake",
  "nakatake\nscrna\nsimulated",
  "joung",
  "norman",
  "replogle1",
  "replogle3",
  "replogle4",
  "adamson",
  "replogle2",
  "replogle2 large effect",
  "replogle2 tf only",
  "replogle2_large_effect",
  "replogle2_tf_only",
  "replogle2\nlarge effect",
  "replogle2\nlarge\neffect",
  "replogle2\ntf only",
  "freimer",
  "dixit", 
  "frangieh_IFNg_v1",
  "frangieh\nIFNg v1",
  "frangieh\nIFNg\nv1",
  "frangieh IFNg v1",
  "frangieh_IFNg_v2",
  "frangieh\nIFNg v2",
  "frangieh\nIFNg\nv2",
  "frangieh IFNg v2",
  "frangieh_IFNg_v3",
  "frangieh\nIFNg v3",
  "frangieh\nIFNg\nv3",
  "frangieh IFNg v3"
)
DEGREE_COLUMNS = c(
  "in-degree_ANANSE_0.5",
  "in-degree_ANANSE_tissue_0.5",
  "in-degree_MARA_FANTOM4",
  "in-degree_STRING",
  "in-degree_cellnet_human_Hg1332",
  "in-degree_cellnet_human_Hugene",
  "in-degree_cellnet_mouse_4302",
  "in-degree_cellnet_mouse_mogene",
  "in-degree_celloracle_human",
  "in-degree_chea",
  "in-degree_csnets",
  "in-degree_encode-nets_human",
  "in-degree_fntm",
  "in-degree_gtex_rna",
  "in-degree_humanbase",
  "in-degree_magnum_compendium_32",
  "in-degree_magnum_compendium_394",
  "in-degree_magnum_compendium_ppi",
  "in-degree_regulatorynetworks.org_human",
  "in-degree_regulatorynetworks.org_mouse",
  "out-degree_ANANSE_0.5",
  "out-degree_ANANSE_tissue_0.5",
  "out-degree_MARA_FANTOM4",
  "out-degree_STRING",
  "out-degree_cellnet_human_Hg1332",
  "out-degree_cellnet_human_Hugene",
  "out-degree_cellnet_mouse_4302",
  "out-degree_cellnet_mouse_mogene",
  "out-degree_celloracle_human",
  "out-degree_chea",
  "out-degree_csnets",
  "out-degree_encode-nets_human",
  "out-degree_fntm",
  "out-degree_gtex_rna",
  "out-degree_humanbase",
  "out-degree_magnum_compendium_32",
  "out-degree_magnum_compendium_394",
  "out-degree_magnum_compendium_ppi",
  "out-degree_regulatorynetworks.org_human",
  "out-degree_regulatorynetworks.org_mouse"
)
EVAL_COLUMNS = c(
  c(
    "condition",
    "group",
    "num_observations_in_group",
    "gene",
    "n_exons",
    "pLI",
    "highly_variable",
    "highly_variable_rank",
    "means",
    "variances",
    "variances_norm",
    "spearman",
    "pearson",
    "mae",
    "mse",
    "mse_top_20",
    "mse_top_100",
    "mse_top_200",
    "overlap_top_20",
    "overlap_top_100",
    "overlap_top_200",
    "pearson_top_20",
    "pearson_top_100",
    "pearson_top_200",
    "proportion_correct_direction",
    "pvalue_effect_direction",
    "pvalue_targets_vs_non_targets",
    "fc_targets_vs_non_targets",
    "expression_level_after_perturbation",
    "perturbation_type",
    "prediction_timescale",
    "is_control",
    "timepoint",
    "cell_type",
    "cell_type_correct",
    "cell_label_accuracy",
    "distance_in_pca",
    "is_timescale_strict",
    "unique_id",
    "factor_varied",
    "regression_method",
    "feature_extraction",
    "predict_self",
    "merge_replicates",
    "perturbation_dataset",
    "num_genes",
    "starting_expression",
    "network_datasets",
    "network_prior",
    "data_split_seed",
    "eligible_regulators",
    "cell_type_sharing_strategy",
    "low_dimensional_structure",
    "low_dimensional_training",
    "low_dimensional_value",
    "matching_method"
  )
)

reorder_datasets = function(datasets){
  datasets = as.character(datasets)
  missing_levels = setdiff(unique(datasets), DATASET_ORDER)
  if(length(missing_levels)>0){
    print("Reordering datasets, found some missing from the hardwired DATASET_ORDER:")
  }
  print(missing_levels)
  datasets %<>% factor(levels = c(DATASET_ORDER, missing_levels))
  return(datasets)
}

#' Load evaluation results from a list of benchmarking experiments.
#'
#' @param experiments
#' @param stratify_by_pert If true (default), collect results at per-perturbation resolution. Otherwise, collect results at per-target-gene resolution.
#'
collect_experiments = function(
    experiments, 
    stratify_by_pert = T
){
  X <- list()
  for (experiment in experiments) {
    if (stratify_by_pert){
      filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerPert.parquet")
    } else {
      filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerTarget.parquet")
    }
   
    try({
      X[[experiment]] <- arrow::read_parquet(filepath, as_data_frame = T, mmap = T)
      if (is.null(X[[experiment]][["cell_type"]])){
        X[[experiment]][["cell_type"]] = 0
      }
      X[[experiment]][["cell_type"]] %<>% as.character
      if (is.null(X[[experiment]][["louvain"]])){
        X[[experiment]][["louvain"]] = 0
      }
      X[[experiment]][["refers_to"]] %<>% as.character
      X[[experiment]][["question"]] %<>% as.character
      X[[experiment]][["louvain"]] %<>% as.character
    })
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
  try({X$perturbation_dataset %<>% gsub("paul.", "paul", .)}, silent = T) # Paul1 and Paul2 are separate for evals but really go together
  try({X$perturbation_dataset %<>% gsub("replogle", "replogle1", .)}, silent = T) # we renamed replogle to replogle1
  try({X$perturbation_dataset %<>% gsub("replogle11", "replogle1", .)}, silent = T) # we renamed replogle to replogle1
  try({X$perturbation_dataset %<>% gsub("replogle12", "replogle2", .)}, silent = T) # we renamed replogle to replogle1
  try({X$perturbation_dataset %<>% gsub("replogle13", "replogle3", .)}, silent = T) # we renamed replogle to replogle1
  try({X$perturbation_dataset %<>% gsub("replogle14", "replogle4", .)}, silent = T) # we renamed replogle to replogle1
  try({X$perturbation_dataset %<>% reorder_datasets}, silent = T)
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

best_method = function(mae, x){
  return(x[which.min(mae)])
}
  
percent_change_from_best_baseline = function(mae, x){
  m = min(mae[x %in% c("mean", "median")], na.rm = T)
  return(100*((m - mae) / m))
}

percent_change_from_best = function(x){
  m = max(x, na.rm = T)
  x = -abs( 100*( (m - x) / (m+0.00000001) ) )
  return(x)
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
    metrics = c(# "pearson_top_20",
                "pearson_top_100", 
                # "pearson_top_200",
                # "overlap_top_20",              
                "overlap_top_100",                 
                # "overlap_top_200",
                # "mse_top_20",
                "mse_top_100", 
                # "mse_top_200",
                "mse", 
                "mae", 
                "spearman", 
                "proportion_correct_direction", 
                "cell_label_accuracy"),
    metrics_where_bigger_is_better = c(
      "spearman", "proportion_correct_direction",                 
      "pearson_top_20", 
      "pearson_top_100", 
      "pearson_top_200", 
      "overlap_top_20",                 
      "overlap_top_100",                 
      "overlap_top_200",                 
      "cell_label_accuracy", 
      "proportion correct direction",
      "cell type correct"), 
    scales = "free", 
    do_wrap = F
){
  X[["facet1"]] = X[[facet1]]
  X[["facet2"]] = X[[facet2]]
  X %<>%
    group_by(x, facet1, facet2) %>%
    summarise(across(metrics, \(v) sum(v*num_observations_in_group)/sum(num_observations_in_group)))
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
    mutate(scaled_value = percent_change_from_best(value)) %>%
    mutate(rank_where_lower_is_better = rank(-scaled_value))
  X %<>% subset(!is.na(scaled_value))
  mean_rank = X %>% 
    group_by(x, facet1, facet2) %>%
    summarize(metric = "MEAN ACROSS ALL METRICS", rank_where_lower_is_better = mean(rank_where_lower_is_better))
  X = rbind(X, mean_rank)
  X[["metric"]] %<>% factor(levels = gsub("_", " ", metrics) %>% c("MEAN ACROSS ALL METRICS")) #order y axis in plot
  # plawt
  g = ggplot(X) +
    geom_tile(aes(x = x, y = metric, fill = rank_where_lower_is_better)) + 
    # scale_fill_gradient( breaks=c(0,1),labels=c("min","max"), limits=c(0,1)) +
    # geom_point(data = subset(X, is_best), aes(x = x, y = metric, color = is_best)) + 
    # scale_color_manual(values = c("red")) +
    # labs(x = "", y = "", fill = "Scaled\nvalue", color = "Best\nperformer") +
    # labs(x = "", y = "", fill = "Percent change\nfrom best \n(capped at 10%)", color = "Best\nperformer") +
    labs(x = "", y = "", fill = "Rank (lower is better)", color = "Best\nperformer") +
    scale_fill_viridis_c(direction = -1) + 
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black")
      ) 
  if (do_wrap){
    g = g + facet_wrap(~facet1, scales = scales) 
  } else {
    g = g + facet_grid(facet1~facet2, scales = scales) 
  }
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


do_subnetworks_match = function(
    perturbation_dataset, 
    network,
    cell_type_matching = rjson::fromJSON(file = "../../accessory_data/cell_type_matching.json")
){
  networks_by_dataset = cell_type_matching$celltypes[perturbation_dataset]
  networks_by_dataset %<>% unlist
  # If x is from a network that does not include any relevant subnetwork, the key might be missing from cell_type_matching$networks.
  # If x is from a network that does include a relevant subnetwork, this code runs as expected: check if this subnetwork is the relevant one.
  try(
    { 
      networks_by_dataset = cell_type_matching$networks[networks_by_dataset][[1]]
      # This incomprehensible one-liner converts a nested list of vectors to a tidy dataframe: {a:[1,2], b:2} becomes [[a,a,b], [1,2,2]].
      networks_by_dataset = data.frame(network = Reduce(c, mapply(rep, names(networks_by_dataset), sapply(networks_by_dataset, length))), 
                     subnetwork = Reduce(c, networks_by_dataset))
      return( network %in% c(networks_by_dataset$network, paste(networks_by_dataset$network, networks_by_dataset$subnetwork)) )
    }, 
    silent = T
  )
  return(F)
}

add_network_cell_type_metadata = function(
    X,   
    single_networks = c("celloracle_human",  
                        "magnum_compendium_ppi" ,
                        "MARA_FANTOM4"     , 
                        "STRING",   
                        "magnum_compendium_32", 
                        "dense", 
                        "empty", 
                        "endoderm" )
){
  X$network_cell_type_matches = mapply(
    do_subnetworks_match, 
    X$perturbation_dataset, 
    X$network_datasets
  )
  X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
  X$perturbation_dataset %<>% gsub("γ", "g", .)
  X %<>% make_the_usual_labels_nice()
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
    gsub("multipotent", "", .) %>%
    gsub("unrestricted", "", .) %>%
    gsub("somatic", "", .) %>%
    gsub("acute myeloid leukemia", "aml", .) %>%
    gsub("peripheral blood mononuclear cells", "pbmc", .) %>%
    gsub("suprapubic", "", .) %>%
    gsub("lower leg", "", .) %>%
    gsub("brain .*", "brain", .) %>%
    gsub("cell line", "", .) %>%
    gsub("muscleskel", "muscle", .) %>%
    gsub("pancreatic", "pancreas", .) 
  X$network_tissue[X$network_source %in% single_networks] = X$network_source[X$network_source %in% single_networks] 
  X$network_source[X$network_source %in% single_networks] = "other"
  X$network_pretty = paste(
    as.integer(as.factor(X$network_source)),
    X$network_tissue
  )
  return(X)
}
