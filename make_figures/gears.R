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
    filepath <- paste0("../experiments/", experiment, "/outputs/evaluationPerPert.parquet")
    X[[experiment]] <- arrow::read_parquet(filepath)
    X[[experiment]]$refers_to %<>% as.character
    X[[experiment]]$question %<>% as.character
  }
  X <- bind_rows(X)
  return(X)
}

X = collect_experiments(c("1.4.2_1","1.4.2_2","1.4.2_3","1.4.2_5","1.4.2_6"))
X %<>%
  group_by(regression_method, eligible_regulators, perturbation_dataset, desired_heldout_fraction) %>%
  summarise(across(ends_with("benefit"), mean))
X$desired_heldout_fraction %<>% 
  multiply_by(100) %>%
  paste0("Held-out : ", ., "%")
for(metric in c("mae_benefit")){
  ggplot(X) + 
    geom_point(aes_string(x = "regression_method", 
                          y = metric,
                          color='eligible_regulators'), position = position_dodge(width=0.3)) + 
    labs(x='', 
         y = "MAE improvement over baseline",
         color='Eligible regulators') +
    facet_grid(desired_heldout_fraction~perturbation_dataset) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}
ggsave(paste0('plots/fig_gears_', metric, '.pdf'), width = 10, height = 5)
