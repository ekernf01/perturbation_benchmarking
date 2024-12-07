library("ggplot2")
library("dplyr")
library("stringr")
library("arrow")
library("magrittr")
library("rjson")
# # Versions of AnnDataR and hdf5r we used:
# remotes::install_version("hdf5r", "1.3.11")
# devtools::install_github("scverse/anndataR", ref = "dbc6897")

setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
source("plotting_functions.R")

# Definitive Endoderm
{
  embeddings = load_embeddings("1.2.2_14")
  # Demo plot of velocity
  ggplot(embeddings %>% subset(perturbation=="SOX17" & matching_method=="optimal_transport")) +
    ggtitle("Lineage relationships predicted by optimal transport") +
    geom_point(aes(x = train_viz1, y = train_viz2),
               color = "gray") + 
    geom_segment(aes(x = progenitor_viz1 , y = progenitor_viz2,
                     xend = train_viz1, yend = train_viz2),
                 color = "blue",
                 alpha = 0.01,
                 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
    theme_minimal() 
  ggsave("timeseries_plots/endoderm_viz_example_velocity.pdf", width = 8, height = 8)
  # Demo plot showing perturbation effects
  ggplot() +
    ggtitle("Predicted perturbation effects (SOX17 knockdown)") +
    geom_point(
      aes(x = train_viz1, y = train_viz2),
      color = "gray", 
      data = embeddings %>% subset( (perturbation == "SOX17") & prediction_timescale==1 ),
      ) + 
    geom_segment(
      aes(x = train_viz1 , y = train_viz2,
          xend = train_viz1 + predicted_delta_viz1, 
          yend = train_viz2 + predicted_delta_viz2),
      color = "red",
      alpha = 0.1,
      data = embeddings %>% subset( abs(predicted_delta_viz1) + abs(predicted_delta_viz2) >= 0.5 ),
      arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
    theme_minimal() + 
    facet_wrap(~matching_method)
  ggsave("timeseries_plots/endoderm_viz_example_SOX17.pdf", width = 8, height = 8)
  # Compare perturbation scores versus CRISPR screen
  overall_scores = embeddings %>%
    group_by(cell_type, timepoint, perturbation, prediction_timescale, matching_method) %>%
    summarize(
      perturbation_score = mean(perturbation_score), 
      li_et_al_biggest_effect = maxabs(c(
        mean(brunello.neg.lfc),
        mean(brunello.pos.lfc),
        mean(gecko.neg.lfc),
        mean(gecko.pos.lfc)
      ))
    ) %>%
    group_by(cell_type, prediction_timescale, matching_method) %>%
    mutate( perturbation_score_rank = rank(-abs(perturbation_score)) ) %>%
    mutate( top_30 = ifelse(perturbation_score_rank<=30, "top 30 predictions", "other") ) %>%
    subset(cell_type != "pluripotent") 
  
  ggplot(overall_scores) + 
    facet_grid(prediction_timescale~cell_type, scale = "free_x") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_point(aes(perturbation_score, li_et_al_biggest_effect, color = matching_method, alpha = top_30)) + 
    ggrepel::geom_label_repel(
      data=subset(overall_scores, perturbation_score_rank<=5) %>% 
        subset(!duplicated(perturbation)), 
      mapping = aes(perturbation_score, li_et_al_biggest_effect, label = perturbation)
    ) +
    ylab("Largest gRNA abundance effect in Li et al.")

  overall_scores %>% 
    subset(perturbation_score_rank<=30) %>% 
    write.csv("timeseries_plots/definitive_endoderm_ps_vs_screen_top_30.csv")
    
  table_to_report %>%
    subset(prediction_timescale==10) %>%
    group_by(cell_type, perturbation, matching_method) %>%
    summarize(
      perturbation_score = maxabs(perturbation_score),
      li_et_al_biggest_effect = maxabs(li_et_al_biggest_effect),
    ) %>% 
    group_by(cell_type, matching_method) %>%
    mutate(is_outlier = 
             rank(-abs(perturbation_score)) <= 30
    ) %>% 
    write.csv("timeseries_plots/endoderm_ps_vs_screen_top_30.csv")
    
  ggplot() + 
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +
    geom_point(aes(x = perturbation_score, y = li_et_al_biggest_effect, color = prediction_timescale), data = overall_scores) + 
    ggrepel::geom_label_repel(aes(x = perturbation_score, y = li_et_al_biggest_effect, label = perturbation), data = outliers) + 
    facet_grid(cell_type~matching_method, scales = "free") + 
    xlab("Cell type-specific perturbation score") +
    ylab("Log fold change in gRNA abundance \n(Largest magnitude from Li et al. 2019 genome-wide screen)") + 
    ggtitle("Predicted and observed effects on endoderm differentiation",
            subtitle = "Labels are shown for the largest absolute values on the x axis.")
  ggsave("timeseries_plots/endoderm_ps_vs_screen.pdf", width = 12, height = 8)
}

{
  embeddings = load_embeddings("1.2.2_15")
  embeddings[embeddings$perturbation=="control","Annotation_summary"] = "Simulated control" 
  ggplot(embeddings %>% subset(perturbation=='control' & prediction_timescale==1 & matching_method=="steady_state")) +
    geom_point(aes(x = train_viz1, y = train_viz2, color = cell_type)) + 
    geom_segment(aes(x = progenitor_viz1 , y = progenitor_viz2,
                     xend = train_viz1, yend = train_viz2),
                 color = "blue",
                 alpha = 0.2,
                 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "open")) +
    theme_minimal() 
  
  embeddings[["super_cell_type"]] = ifelse( 
    embeddings[["cell_type"]] %in% c(
      "Erythroids",
      "Megakaryocytes",
      "MEP"    
    ),
    "ME",
    "GM"
  )
  overall_scores = embeddings %>% 
    group_by(super_cell_type, perturbation, prediction_timescale, matching_method) %>%
    summarize(
      perturbation_score = mean(perturbation_score, na.rm = T), 
      lit_review = unique(Annotation_summary)
    ) %>%   
    ungroup() %>%
    tidyr::pivot_wider(names_from = "super_cell_type", values_from = "perturbation_score") 
  
  outliers = overall_scores %>% 
    group_by(matching_method) %>%
    mutate(is_outlier = rank(-abs(GM)) <= 10 | rank(-abs(ME)) <= 10 ) %>% 
    subset(is_outlier)
  
  overall_scores %>% 
    ggplot() + 
    geom_point(aes(x = GM, y = ME, color = lit_review)) +
    scale_color_manual(values = c(
      "GM"="green", 
      "ME"="red", 
      "GM & ME, stemness"="yellow", 
      "No known role"= "black", 
      "Not reviewed" = "gray", 
      "Simulated control" = "blue")
      ) +
    ggrepel::geom_text_repel(aes(x = GM, y = ME, label = perturbation), data = outliers) + 
    facet_wrap(~matching_method, scales = "free")
  ggsave("timeseries_plots/mouse_blood_ps_vs_screen.pdf", width = 12, height = 8)
  
  table_to_report = overall_scores %>% 
    group_by(matching_method) %>%
    mutate(is_outlier = rank(-abs(GM)) <= 30 | rank(-abs(ME)) <= 30 ) %>% 
    subset(is_outlier) 
  table_to_report %>% subset(is_outlier, select = c("perturbation", "lit_review")) %>% distinct() %>% extract2(2) %>% table
  table_to_report %>% 
    write.csv("timeseries_plots/mouse_blood_ps_vs_screen_top_30.csv")
}

{
  embeddings = load_embeddings("1.2.2_18")
  ggplot(embeddings %>% subset(perturbation=='control' & prediction_timescale==1)) +
    geom_point(aes(x = train_viz1, y = train_viz2, color = cell_type)) + 
    geom_segment(aes(x = progenitor_viz1 , y = progenitor_viz2,
                     xend = train_viz1, yend = train_viz2),
                 color = "blue",
                 alpha = 0.05,
                 arrow = arrow(length=unit(0.30,"cm"), ends="last", type = "closed")) +
    theme_minimal()
  overall_scores = embeddings %>% 
    group_by(cell_type, timepoint, perturbation, prediction_timescale, matching_method) %>%
    summarize(
      perturbation_score = mean(perturbation_score), 
      lit_review = unique(Reference..Pubmed.ID.)
    ) %>%
    tidyr::pivot_wider(names_from = "cell_type", values_from = "perturbation_score")
  ggplot(overall_scores) + 
    geom_point(aes(x = notochord, y = `mesodermal progenitor cells (contains PSM)`, color = timepoint))
  ggsave("timeseries_plots/axial_mesoderm_ps_vs_screen.pdf", width = 12, height = 8)
  table_to_report = overall_scores %>% 
    group_by(matching_method) %>%
    mutate(is_outlier = rank(-abs(`mesodermal progenitor cells (contains PSM)`)) <= 30 | rank(-abs(notochord)) <= 30 ) %>% 
    subset(is_outlier) 
  table_to_report %>% subset(is_outlier, select = c("perturbation", "lit_review")) %>% distinct() %>% table
  table_to_report %>% 
    write.csv("timeseries_plots/mouse_blood_ps_vs_screen_top_30.csv")
}

# Numbers 16-21 are other datasets that are simpler to deal with. 
{
  X = collect_experiments(paste0("1.2.2_", 16:21)) %>% make_the_usual_labels_nice
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, matching_method) %>%
    summarise(across(DEFAULT_METRICS, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(DEFAULT_METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = DEFAULT_METRICS)
  X[["prediction_timescale"]] %<>% as.character()
  X[["prediction_timescale"]] %<>% factor(levels = gtools::mixedsort(unique(X[["prediction_timescale"]])))
  plot = ggplot(X) +
    geom_bar(aes(x = prediction_timescale, y = value, fill = matching_method), position = "dodge", stat = "identity") +
    facet_grid(metric~perturbation_dataset, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
    ylab("")
  print(plot)
  ggsave('timeseries_plots/matching.svg', plot, width = 4, height = 8)
}

{
  X = collect_experiments(paste0("1.4.1_", 0:4))
  X %<>% subset(is_timescale_strict)
  X[X[["pearson_top_20"]]!=0, ]
  X %<>% add_network_cell_type_metadata
  X %<>% make_the_usual_labels_nice
  X %<>% subset(is_timescale_strict)
  X %<>% subset(prediction_timescale==10)
  X$pearson %<>% tidyr::replace_na(0)
  X$pearson_top_20 %<>% tidyr::replace_na(0)
  X$pearson_top_100 %<>% tidyr::replace_na(0)
  X$pearson_top_200 %<>% tidyr::replace_na(0)
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, network_datasets, network_source, network_tissue, network_pretty, network_cell_type_matches ) %>%
    summarise(across(DEFAULT_METRICS, mean, na.rm = T))
  X %<>% tidyr::pivot_longer(cols = all_of(DEFAULT_METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gtools::mixedsort(unique(X[["metric"]])))
  X %<>% subset(network_source != "humanbase")
  X$source = X$network_source %>% 
    gsub("_tissue", "", . ) %>% 
    gsub("cellnet_human.*", "cellnet_human", .) %>% 
    gsub("cellnet_mouse.*", "cellnet_mouse", .)
  for(mymetric in c("overlap top 20", "pearson top 20")){
    plot = X %>% 
      subset(metric == mymetric) %>%
      ggplot() + 
      geom_bar(
        aes(
          x = network_pretty, 
          y = value, 
          color = source, 
          fill = c("Different", "Matched")[network_cell_type_matches+1]
        ), 
        stat = "identity", 
        position = "dodge") +
      scale_fill_manual(values = c("gray", "black")) +
      labs(fill = "Network cell type matches \nperturbation data cell type?", color = "Network source") +
      facet_wrap(~perturbation_dataset, scales = "free", ncol = 1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
      xlab("") + 
      ylab(mymetric)
    print(plot)
    ggsave(paste0('timeseries_plots/cell_type_specific_', mymetric, '.svg'), plot, width = 12, height = 12)
  }
}



