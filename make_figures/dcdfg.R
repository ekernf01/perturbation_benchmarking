library(ggplot2)
library(dplyr)
library(arrow)
library(magrittr)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")

collect_experiments = function(experiments){
  X <- list()
  for (experiment in experiments) {
    filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerPert.parquet")
    X[[experiment]] <- arrow::read_parquet(filepath)
    X[[experiment]]$refers_to %<>% as.character
    X[[experiment]]$question %<>% as.character
  }
  X <- bind_rows(X)
  return(X)
}
X = collect_experiments(c("1.6.1_1", "1.6.1_3", "1.6.1_6", "1.6.1_7","1.6.1_8", "1.6.1_9", "1.6.1_10", "1.6.1_11"))
X <- X %>% mutate(chart_x = paste(regression_method, starting_expression, sep = "_"))
method_tidy = c(
  "DCDFG-spectral_radius-mlplr-False"="DCD-FG" ,
  "median"  = "median",                              
  "DCDFG-spectral_radius-linearlr-False"="NOTEARS-LR",
  "mean"  = "mean"   
) 
X$regression_method = method_tidy[X$regression_method] %>% factor(levels = c("median", "mean", "NOTEARS-LR", "DCD-FG"))
X$perturbation_dataset %<>% gsub("γ", "g", .)
X %<>%
  group_by(regression_method, starting_expression, perturbation_dataset) %>%
  summarise(across(ends_with("benefit"), mean))
my_levels = unique(X$perturbation_dataset)
X$perturbation_dataset %<>% factor(levels = unique(c("frangieh_IFNg_v1", "frangieh_IFNg_v3", "nakatake", my_levels)))
for (metric in c("mae_benefit")) {
  ggplot(X, aes(x = regression_method, y = mae_benefit, fill = starting_expression, color = starting_expression)) +
    geom_point(position = position_dodge(width = 0.5)) +
    facet_wrap(~perturbation_dataset, scales = "free", ncol = 4) +
    labs(x = "Method",
         y = "MAE improvement over baseline") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  
  dir.create("plots", showWarnings = FALSE)
  ggsave(filename = paste0("plots/fig_dcdfg_", metric, ".pdf"), width = 8, height = 4)
}

