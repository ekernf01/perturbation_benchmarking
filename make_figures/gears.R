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

X = collect_experiments(c(
  "1.4.2_1",
  "1.4.2_2",
  "1.4.2_3",
  # "1.4.2_4",
  "1.4.2_5",
  "1.4.2_6",
  "1.4.2_7"
  # "1.4.2_8"
)
)
X$regression_method %<>% gsub("0$", "", .)
the_usual_levels = unique(X$regression_method)
X$regression_method %<>% factor(levels = unique(c("empty", "dense", "median", "mean", "celloracle human", the_usual_levels)))

X %<>%
  group_by(regression_method, eligible_regulators, perturbation_dataset, desired_heldout_fraction, data_split_seed) %>%
  summarise(across(c("mae", "mse_top_20"), mean))
X$desired_heldout_fraction %<>% 
  multiply_by(100) %>%
  paste0("Held-out : ", ., "%")
X$data_split_seed %<>% as.character()
for(metric in c("mae", "mse_top_20")){
  ggplot(X) + 
    geom_line(aes_string(x = "regression_method", 
                        y = metric, 
                        group = "data_split_seed",
                        color = "data_split_seed")) + 
    labs(x='', 
         y = metric) +
    facet_grid(perturbation_dataset~desired_heldout_fraction, scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
  ggsave(paste0('plots/fig_gears_', metric, '.pdf'), width = 8, height = 3)
}


