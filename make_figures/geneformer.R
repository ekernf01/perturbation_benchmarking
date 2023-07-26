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

X = collect_experiments(c("1.3.3_1", "1.3.3_2", "1.3.3_3", "1.3.3_5", "1.3.3_6", "1.3.3_7", "1.3.3_8", "1.3.3_9", "1.3.3_10"))
X$regression_method %<>% gsub("0$", "", .)
X %<>%
  group_by(regression_method, eligible_regulators, perturbation_dataset, desired_heldout_fraction) %>%
  summarise(across(ends_with("benefit"), mean))
X$desired_heldout_fraction %<>% 
  multiply_by(100) %>%
  paste0("Held-out : ", ., "%")
for(metric in c("mae")){
  ggplot(X) + 
    geom_point(aes_string(x = "regression_method", 
                          y = metric), position = position_dodge(width=0.3)) + 
    labs(x='', 
         y = "Mean absolute error") +
    facet_wrap(~perturbation_dataset, scales = "free_y", nrow = 2) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
}
ggsave(paste0('plots/fig_geneformer_', metric, '.pdf'), width = 8, height = 3)
