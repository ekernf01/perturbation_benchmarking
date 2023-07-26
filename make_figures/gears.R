library(ggplot2)
library(dplyr)
library(stringr)
library(arrow)
library(magrittr)
setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")

# Figure <gears>
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

X = collect_experiments(c("1.4.2_1","1.4.2_2","1.4.2_3","1.4.2_5","1.4.2_6","1.4.2_7"))
X$regression_method %<>% gsub("0$", "", .)
the_usual_levels = unique(X$regression_method)
X$regression_method %<>% factor(levels = unique(c("empty", "dense", "median", "mean", "celloracle human", the_usual_levels)))

X %<>%
  group_by(regression_method, eligible_regulators, perturbation_dataset, desired_heldout_fraction) %>%
  summarise(across(starts_with("mae"), mean))
X$desired_heldout_fraction %<>% 
  multiply_by(100) %>%
  paste0("Held-out : ", ., "%")
for(metric in c("mae")){
  ggplot(X) + 
    geom_point(aes_string(x = "regression_method", 
                          y = metric,
                          color='eligible_regulators'), position = position_dodge(width=0.3)) + 
    labs(x='', 
         y = "Mean absolute error",
         color='Eligible regulators') +
    facet_grid(desired_heldout_fraction~perturbation_dataset) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}
ggsave(paste0('plots/fig_gears_', metric, '.pdf'), width = 8, height = 3)


