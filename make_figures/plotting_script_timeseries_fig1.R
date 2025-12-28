# environment setup
{
  library("viridis")
  library("gtools")
  library("ggplot2")
  library("dplyr")
  library("stringr")
  library("arrow")
  library("magrittr")
  library("rjson")
  library("ggrepel")
  setwd("/home/ekernf01/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/perturbation_benchmarking/make_figures/")
  source("plotting_functions.R")
  METRICS = "pearson"
}

# ===================== Figure 1, perturbation scores and lit reviews ====================================
# timescale experiments
{
  X = collect_experiments(experiments = paste0("1.2.2_", c(14:17))) %>% make_the_usual_labels_nice
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, matching_method) %>%
    summarise(across(METRICS, mean, na.rm = T))
  X %<>% tidyr::pivot_longer(cols = all_of(METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gsub("_", " ", METRICS))
  X[["prediction_timescale"]] %<>% as.character()
  X[["prediction_timescale"]] %<>% factor(levels = gtools::mixedsort(unique(X[["prediction_timescale"]])))
  plot = ggplot(X) + 
    geom_bar(aes(x = prediction_timescale, y = value, fill = matching_method), position = "dodge", stat = "identity") +
    facet_grid(1~perturbation_dataset, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
    ylab("Pearson correlation between predicted\nand observed differential expression")
  print(plot)
  ggsave('timeseries_plots/matching.svg', plot, width = 8, height = 3)
}
# Definitive Endoderm
{
  embeddings = load_embeddings("1.2.2_14")
  lit_review = read.csv("../../perturbation_data/perturbations/definitive_endoderm/lit_review.csv")
  lit_review %<>% rename(cell_type = Cell.type.affected)
  lit_review %<>% mutate(cell_type = gsub("primitive streak", "mesendoderm", tolower(cell_type)))
  lit_review %<>% tidyr::complete(perturbation, cell_type, fill = list(effect_direction="none"))
  embeddings %<>% merge(lit_review, by = c("perturbation", "cell_type"), all.x = T)
  embeddings %<>% mutate(Included.in.literature.review. = ifelse(is.na(Included.in.literature.review.), "no", Included.in.literature.review.))
  embeddings %<>% mutate(effect_direction = ifelse(is.na(effect_direction), "none", effect_direction))
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
  ggsave("timeseries_plots/definitive_endoderm_viz_example_velocity.pdf", width = 8, height = 8)
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
  ggsave("timeseries_plots/definitive_endoderm_viz_example_SOX17.pdf", width = 8, height = 8)
  # Compare perturbation scores versus CRISPR screen
  overall_scores = embeddings %>%
    group_by(cell_type, perturbation, prediction_timescale, matching_method) %>%
    summarize(
      perturbation_score = mean(perturbation_score), 
      li_et_al_biggest_effect = maxabs(c(
        mean(brunello.neg.lfc),
        mean(brunello.pos.lfc),
        mean(gecko.neg.lfc),
        mean(gecko.pos.lfc)
      )),
      Included.in.literature.review. = unique(Included.in.literature.review.), 
      effect_direction = unique(effect_direction),
      PMID = unique(PMID),
      Notes = unique(Notes)
    ) %>%
    group_by(cell_type, prediction_timescale, matching_method) %>%
    mutate( perturbation_score_rank = rank(-abs(perturbation_score)) ) %>%
    mutate( top_30 = ifelse(perturbation_score_rank<=30, "top 30 predictions", "other") ) %>%
    subset(cell_type != "pluripotent") 
  # perturbation scores versus CRISPR screen
  ggplot(overall_scores) + 
    facet_grid(~cell_type, scale = "free_x") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_point(aes(perturbation_score, li_et_al_biggest_effect, color = matching_method, alpha = top_30)) + 
    ggrepel::geom_text_repel(
      data=subset(overall_scores, perturbation_score_rank<=5) %>% 
        group_by(perturbation, cell_type) %>%
        summarise(perturbation_score = maxabs(perturbation_score), 
                  li_et_al_biggest_effect = li_et_al_biggest_effect[[1]]),
      mapping = aes(perturbation_score, li_et_al_biggest_effect, label = perturbation)
    ) +
    xlab("Cell type-specific perturbation score") +
    ylab("Log fold change in gRNA abundance \n(Largest magnitude from Li et al. 2019 genome-wide screen)") + 
    ggtitle("Predicted and observed effects on endoderm differentiation",
            subtitle = "Labels are shown for the largest absolute values on the x axis.")
  ggsave("timeseries_plots/definitive_endoderm_ps_vs_screen.pdf", width = 8, height = 3)
  
  # perturbation scores overall distribution
  overall_scores %>% 
    ggplot() + 
    stat_ecdf(aes(x = perturbation_score)) + 
    facet_grid(matching_method~cell_type) + 
    labs(fill = "Exact zero", y = "Cumulative density") + 
    geom_vline(xintercept = 0, color = "red") + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave("timeseries_plots/definitive_endoderm_ps_distribution.pdf", width = 3, height = 6)
  
  
  # Perturbation scores versus each other
  prediction_discordance = overall_scores %>% 
    subset(perturbation_score_rank<=30) %>% 
    mutate(predicted_effect = sign(perturbation_score)) %>% 
    mutate(method = matching_method %>% gsub("_", " ", .)) %>% 
    extract(c("method", "cell_type", "perturbation", "predicted_effect")) %>% 
    ungroup() %>% 
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect=0)) %>%
    tidyr::pivot_wider(names_from = c("method"), id_cols = c("perturbation", "cell_type"), values_from = "predicted_effect")
  # Get the correlation matrix and unroll it, keeping rownames. Sorry for the ugly-ass code.
  get_corr = function(current_cell_type) {
    prediction_discordance[prediction_discordance$cell_type==current_cell_type, -(1:2)] %>%
      cor %>%
      as.data.frame %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(cols = -1) %>%
      cbind(cell_type=current_cell_type ) %>%
      set_colnames(c("method1", "method2", "correlation", "cell_type"))
  }
  rbind(get_corr("mesendoderm"), get_corr("endoderm")) %>%
  ggplot() + 
    geom_tile(aes(x = method1, y = method2, fill = correlation)) + 
    facet_wrap(~cell_type, ncol=1) + 
    labs(fill = "Correlation between\npredictions") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0
    ) + 
    coord_fixed() + 
    ggtitle("Comparison of top predictions") + 
    xlab("") + ylab("")
  ggsave("timeseries_plots/definitive_endoderm_discordance_top_30.pdf", width = 4, height = 6) 
  
  
  # perturbation scores versus lit review
  heatmap_data %>% 
    write.csv("timeseries_plots/definitive_endoderm_ps_vs_screen_top_30.csv")
  heatmap_data %>% 
    subset(T, select = c("cell_type", "effect_direction")) %>%
    table() %>%
    write.csv("timeseries_plots/definitive_endoderm_ps_vs_screen_top_30_unique.csv")
  heatmap_data %>% 
    mutate(prediction_for = paste0("Prediction about:\n", cell_type)) %>%
    subset(predicted_effect!="none") %>%
    mutate(accuracy = "Incorrect (no effect)") %>%
    mutate(accuracy = ifelse(effect_direction != "none",           "Incorrect (wrong direction)", accuracy)) %>%
    mutate(accuracy = ifelse(effect_direction == predicted_effect, "Correct",         accuracy)) %>%
    mutate(accuracy = forcats::fct_relevel(accuracy, rev(sort(unique(accuracy)))[c(1,3,2)])) %>%
    ggplot() + 
    geom_bar(aes(x=method, fill = accuracy)) + 
    scale_fill_manual(values = c( "Correct"="green", "Incorrect (wrong direction)" = "red", "Incorrect (no effect)"="gray")) +
    facet_grid(prediction_for~"") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Number of genes") +
    ggtitle("Literature review of top 30 \npredictions for each model")
  ggsave("timeseries_plots/definitive_endoderm_lit_review_top_30.pdf", width = 4, height = 6)
}

# Paul et al. 2015 myelopoiesis data
{
  embeddings = load_embeddings("1.2.2_15")
  embeddings[embeddings$perturbation=="control","Annotation_summary"] = "Simulated control" 
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
      lit_review = unique(Annotation_summary), 
      included_co = unique(Included.in.CellOracle.lit.review.)
    ) %>%   
    ungroup() %>%
    tidyr::pivot_wider(names_from = "super_cell_type", values_from = "perturbation_score") %>% 
    group_by(matching_method, prediction_timescale) %>%
    mutate(GM_abs_rank = rank(-abs(GM)), ME_abs_rank = rank(-abs(ME)) ) 
  prediction_discordance = overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>%
    mutate(GM = sign(GM)) %>%
    mutate(ME = sign(ME)) %>%
    tidyr::pivot_longer(cols = c("GM", "ME"), names_to = "cell_type", values_to = "predicted_effect") %>%
    mutate(method = matching_method %>% gsub("_", " ", .)) %>%
    extract(c("method", "cell_type", "perturbation", "predicted_effect")) %>% 
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect=0)) %>%
    tidyr::pivot_wider(names_from = c("method"), id_cols = c("perturbation", "cell_type"), values_from = "predicted_effect")
  
  get_corr = function(current_cell_type) {
    prediction_discordance[prediction_discordance$cell_type==current_cell_type, -(1:2)] %>%
      cor %>%
      as.data.frame %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(cols = -1) %>%
      cbind(cell_type=current_cell_type ) %>%
      set_colnames(c("method1", "method2", "correlation", "cell_type"))
  }
  # Direction of predicted effect
  rbind(get_corr("GM"), get_corr("ME")) %>%
  ggplot() + 
    geom_tile(aes(x = method1, y = method2, fill = correlation)) + 
    facet_wrap(~cell_type, ncol=1) + 
    labs(fill = "Correlation between\npredictions") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0
    ) + 
    coord_fixed() + 
    ggtitle("Comparison of top predictions") + 
    xlab("") + ylab("")
  ggsave("timeseries_plots/paul_discordance_top_30.pdf", width = 4, height = 7)
  
  # Predictions vs lit review
  overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>% 
    write.csv("timeseries_plots/paul_ps_vs_screen_top_30.csv")
  overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30, select = c("perturbation", "lit_review", "included_co")) %>%
    distinct() %>%
    extract(c("lit_review", "included_co")) %>%
    table() %>%
    write.csv("timeseries_plots/paul_ps_vs_screen_top_30_unique.csv")
  rbind(get_corr("mesendoderm"), get_corr("endoderm"))
  
  
  overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>% 
    subset(abs(GM)>0|abs(ME)>0) %>% 
    mutate(predicted_effect_on_GM = GM_abs_rank<=30, predicted_effect_on_ME = ME_abs_rank<=30) %>%
    tidyr::pivot_longer( cols = c("predicted_effect_on_GM", "predicted_effect_on_ME"), names_to = "cell_type", values_to = "predicted_effect" ) %>% 
    subset(predicted_effect) %>%
    mutate(cell_type = gsub("predicted_effect_on_", "", cell_type)) %>%
    mutate(cell_type = paste0("Prediction about:\n", cell_type)) %>%
    extract(c("matching_method", "lit_review", "prediction_timescale", "cell_type")) %>%
    rename(known_role = lit_review) %>%
    mutate(known_role = ifelse(is.na(known_role), "Not reviewed", known_role)) %>%
    mutate(known_role = forcats::fct_relevel(known_role, "Not reviewed", after = 0)) %>%
    mutate(known_role = forcats::fct_relevel(known_role, "No known role", after = 0)) %>%
    table() %>%
    as.data.frame() %>%
    ggplot() + 
    geom_bar(aes(x=matching_method,y=Freq, fill = known_role), stat="identity") + 
    scale_fill_manual(values = c(
      "GM"="green",
      "ME"="red",
      "GM & ME, stemness"="yellow",
      "DC" = "blue",
      "No known role"= "black",
      "Not reviewed" = "gray")
      ) +
    facet_grid(cell_type~"", scale = "free_y") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Number of genes") +
    labs(fill="Known roles") +
    ggtitle("Literature review of top 30 \npredictions for each model")
  ggsave("timeseries_plots/paul_lit_review_top_30.pdf", width = 4, height = 6)
}

{
  # fantom4 thp-1 data
  embeddings = load_embeddings("1.2.2_16")
  # Compare perturbation scores versus CRISPR screen
  overall_scores = embeddings %>%
    subset(timepoint!=0) %>%
    group_by(cell_type, timepoint, perturbation, prediction_timescale, matching_method) %>%
    summarize(
      perturbation_score = mean(perturbation_score), 
      gRNA_log2_fc = mean(Log2.fold.change),                                       
      FDR.P.value = mean(FDR.P.value),                                       
      Cell.type.affected = unique(Cell.type.affected),
      Known.regulation.of.macrophage.differentiation.or.CD14 = unique(Known.regulation.of.macrophage.differentiation.or.CD14),
      screened_by_surdziel = unique(screened_by_surdziel)
    ) %>%
    group_by(timepoint, prediction_timescale, matching_method) %>%
    mutate( Cell.type.affected = ifelse(Cell.type.affected=="", "none", Cell.type.affected))  %>%
    mutate( perturbation_score_rank = rank(-abs(perturbation_score)) ) %>%
    mutate( top_30 = ifelse(perturbation_score_rank<=30, "top 30 predictions", "other") ) 
  
  # Plot top predictions versus screens
  overall_scores %>% 
    subset(perturbation_score_rank<=30) %>% 
    write.csv("timeseries_plots/fantom4_ps_vs_screen_top_30.csv")
  overall_scores %>% 
    subset(perturbation_score_rank<=30, select = c("perturbation", "Cell.type.affected", "screened_by_surdziel")) %>%
    distinct() %>%
    extract(c("Cell.type.affected", "screened_by_surdziel")) %>%
    table() %>%
    write.csv("timeseries_plots/fantom4_ps_vs_screen_top_30_unique.csv")
  overall_scores %>%
    subset(perturbation_score_rank<=30) %>% 
    mutate(screened_by_surdziel = paste0(screened_by_surdziel, " by \nSurdziel et al.")) %>%
    ggplot() + 
    geom_tile(aes(x = matching_method, 
                  y = reorder(perturbation, -rank(paste0(Cell.type.affected, -rank(screened_by_surdziel)))), 
                  fill = Cell.type.affected)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(fill = "Cell type affected\nin CD14 screen") +
    ylab("Perturbed gene") + 
    xlab("") + 
    facet_wrap(~screened_by_surdziel, scale = "free_y") 
  ggsave("timeseries_plots/thp1_screen_vs_top_30.pdf", width = 5, height = 8)
}

