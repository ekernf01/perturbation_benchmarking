# environment setup
{
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
  METRICS = c(METRICS, "distance_in_pca")
}

# ===================== Figure 1, perturbation scores and lit reviews ====================================
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
  ggsave("timeseries_plots/definitive_endoderm_ps_vs_screen.pdf", width = 8, height = 6)
  
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
  prettify = function(x) x %>% gsub("-1", "LoF inhibits differentiation", .) %>% gsub("1", "LoF promotes differentiation", .) %>% gsub("0", "none", .)
  # Ordering heatmaps is always a nightmare.
  # This makes a rough ordering according to th lit review results. 
  heatmap_data = overall_scores %>% 
    subset(perturbation_score_rank<=30) %>% 
    mutate(predicted_effect = prettify(sign(perturbation_score))) %>%
    mutate(method = matching_method %>% gsub("_", " ", .)) %>%
    extract(c("method", "cell_type", "perturbation", "predicted_effect")) %>% 
    ungroup() %>%
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect="none")) %>%
    merge(lit_review, by = c("cell_type", "perturbation"), all.x = T, all.y = F)
  genes_in_order = heatmap_data %>% 
    dplyr::arrange(effect_direction) %>%
    extract2("perturbation") %>%
    unique()
  heatmap_data %<>% 
    mutate(perturbation=factor(perturbation, levels = genes_in_order)) %>%
    mutate(predicted_effect = gsub("none", "Not in top 30", predicted_effect))
  ggplot(heatmap_data) + 
    geom_tile(aes(x = method, y = perturbation, fill = predicted_effect)) + 
    geom_tile(aes(x = " Literature review", y = perturbation, fill = effect_direction)) + 
    scale_fill_manual(
      values = c(
        "none" = "gray30",
        "Not in top 30" = "gray", 
        "GoF promotes differentiation" = "blue", 
        "LoF inhibits differentiation" = "blue",
        "LoF affects expression state" = "cyan",
        "LoF promotes differentiation" = "yellow"
      ), 
      na.value = "black"
    ) +
    facet_wrap(~cell_type) + 
    labs(fill = "Direction of effect") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    coord_fixed() + 
    ggtitle("Literature review of top 30 \npredictions for each model")
  ggsave("timeseries_plots/definitive_endoderm_discordance_top_30.pdf", width = 5, height = 7)
  
  
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
  
  # Direction of predicted effect
  prettify = function(x) x %>% gsub("-1", "-", .) %>% gsub("1", "+", .) %>% gsub("0", "none", .)
  genes_in_order = overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>%
    group_by(perturbation) %>%
    summarise(rank_by = abs(sum(sign(GM)))+abs(sum(sign(ME)))) %>%
    dplyr::arrange(rank_by) %>%
    extract2("perturbation")
  overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>%
    mutate(perturbation=factor(perturbation, levels = genes_in_order)) %>%
    mutate(GM = prettify(sign(GM))) %>%
    mutate(ME = prettify(sign(ME))) %>%
    tidyr::pivot_longer(cols = c("GM", "ME"), names_to = "cell_type", values_to = "predicted_effect") %>%
    mutate(method = matching_method %>% gsub("_", " ", .)) %>%
    extract(c("method", "cell_type", "perturbation", "predicted_effect")) %>% 
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect="none")) %>%
    ggplot() + 
    geom_tile(aes(x = method, y = perturbation, fill = predicted_effect)) + 
    scale_fill_manual(values = c("-"="blue", "+"="yellow", "none" = "gray"), na.value = "gray") +
    facet_wrap(~cell_type) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
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

{
  # Saunders zebrafish axial mesoderm -- currently not used in figure
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

# ======================== Figures 2 and 4: evaluation of log FC =====================================
# timescale experiments
{
  METRICS = c(METRICS, "distance_in_pca")
  X = collect_experiments(paste0("1.2.2_", 16:21)) %>% make_the_usual_labels_nice
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, matching_method) %>%
    summarise(across(METRICS, mean))
  X %<>% tidyr::pivot_longer(cols = all_of(METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gsub("_", " ", METRICS))
  X[["prediction_timescale"]] %<>% as.character()
  X[["prediction_timescale"]] %<>% factor(levels = gtools::mixedsort(unique(X[["prediction_timescale"]])))
  plot = ggplot(X) +
    geom_bar(aes(x = prediction_timescale, y = value, fill = matching_method), position = "dodge", stat = "identity") +
    facet_grid(metric~perturbation_dataset, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
    ylab("")
  print(plot)
  ggsave('timeseries_plots/matching.svg', plot, width = 8, height = 8)
}
# published methods
dir.create("timeseries_plots/published_methods")
{
  X = collect_experiments(paste0("1.5.1_", 0:2)) %>% make_the_usual_labels_nice
  X %<>%
    group_by(regression_method, perturbation_dataset, prediction_timescale) %>%
    summarise(across(METRICS, mean, na.rm = T))
  X %<>% tidyr::pivot_longer(cols = all_of(METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gsub("_", " ", METRICS))
  plot = ggplot(X) +
    geom_bar(aes(x = regression_method, y = value, fill = factor(as.character(prediction_timescale), levels = c(0:5, 10)), group = prediction_timescale), position = "dodge", stat = "identity") +
    facet_grid(metric~perturbation_dataset, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
    ylab("") + 
    geom_hline(yintercept=0) + 
    labs(fill="Simulation\nduration") + 
    scale_fill_manual(values = viridis::viridis(6), guide = "legend")
  print(plot)
  ggsave('timeseries_plots/published_methods/matching.svg', plot, width = 8, height = 8)
}


# =========================== Figure 3, perturbation scores and lit reviews for published methods ==========
# Definitive Endoderm
{
  embeddings = load_embeddings("1.5.1_0") %>% make_the_usual_labels_nice()
  lit_review = read.csv("../../perturbation_data/perturbations/definitive_endoderm/lit_review.csv")
  lit_review %<>% rename(cell_type = Cell.type.affected)
  lit_review %<>% mutate(cell_type = gsub("primitive streak", "mesendoderm", tolower(cell_type)))
  lit_review %<>% tidyr::complete(perturbation, cell_type, fill = list(effect_direction="none"))
  embeddings %<>% merge(lit_review, by = c("perturbation", "cell_type"), all.x = T)
  embeddings %<>% mutate(Included.in.literature.review. = ifelse(is.na(Included.in.literature.review.), "no", Included.in.literature.review.))
  embeddings %<>% mutate(effect_direction = ifelse(is.na(effect_direction), "none", effect_direction))
  
  # Compare perturbation scores versus CRISPR screen
  overall_scores = embeddings %>%
    group_by(cell_type, perturbation, prediction_timescale, regression_method) %>%
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
    group_by(cell_type, prediction_timescale, regression_method) %>%
    mutate( perturbation_score_rank = rank(-abs(perturbation_score)) ) %>%
    mutate( top_30 = ifelse(perturbation_score_rank<=30, "top 30 predictions", "other") ) %>%
    subset(cell_type != "pluripotent") 
  
  # perturbation scores versus CRISPR screen
  ggplot(overall_scores) + 
    facet_grid(~cell_type, scale = "free_x") +
    geom_hline(yintercept=0) +
    geom_vline(xintercept=0) +
    geom_point(aes(perturbation_score, li_et_al_biggest_effect, color = regression_method, alpha = top_30)) + 
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
  ggsave("timeseries_plots/published_methods/definitive_endoderm_ps_vs_screen.pdf", width = 8, height = 6)
  
  # perturbation scores overall distribution
  overall_scores %>% 
    ggplot() + 
    stat_ecdf(aes(x = perturbation_score, color = factor(as.character(prediction_timescale), levels = c(1:5, 10)), group = prediction_timescale)) + 
    facet_grid(regression_method~cell_type) + 
    labs(y = "Cumulative density") + 
    geom_vline(xintercept = 0, color = "red") + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    labs(color="Simulation\nduration") + 
    scale_color_manual(values = viridis::viridis(6), guide = "legend")
  ggsave("timeseries_plots/published_methods/definitive_endoderm_ps_distribution.pdf", width = 3, height = 6)
  
  
  # Perturbation scores versus each other
  prettify = function(x) x %>% gsub("-1", "LoF inhibits differentiation", .) %>% gsub("1", "LoF promotes differentiation", .) %>% gsub("0", "none", .)
  # Ordering heatmaps is always a nightmare.
  # This makes a rough ordering according to th lit review results. 
  heatmap_data = overall_scores %>% 
    subset(perturbation_score_rank<=30) %>% 
    mutate(predicted_effect = prettify(sign(perturbation_score))) %>%
    mutate(method = regression_method %>% gsub("_", " ", .)) %>%
    mutate(method = interaction(regression_method, prediction_timescale)) %>%
    extract(c("method", "regression_method", "prediction_timescale", "cell_type", "perturbation", "predicted_effect", "prediction_timescale")) %>% 
    ungroup() %>%
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect="none")) %>%
    merge(lit_review, by = c("cell_type", "perturbation"), all.x = T, all.y = F)
  genes_in_order = heatmap_data %>% 
    dplyr::arrange(effect_direction) %>%
    extract2("perturbation") %>%
    unique()
  heatmap_data %<>% 
    mutate(perturbation=factor(perturbation, levels = genes_in_order)) %>%
    mutate(predicted_effect = gsub("none", "Not in top 30", predicted_effect))
  ggplot(heatmap_data) + 
    geom_tile(aes(x = method, y = perturbation, fill = predicted_effect)) + 
    geom_tile(aes(x = " Literature review", y = perturbation, fill = effect_direction)) + 
    scale_fill_manual(
      values = c(
        "none" = "gray30",
        "Not in top 30" = "gray", 
        "GoF promotes differentiation" = "blue", 
        "LoF inhibits differentiation" = "blue",
        "LoF affects expression state" = "cyan",
        "LoF promotes differentiation" = "yellow"
      ), 
      na.value = "black"
    ) +
    facet_wrap(~cell_type) + 
    labs(fill = "Direction of effect", x = "Method and prediction timescale") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    coord_fixed() + 
    ggtitle("Literature review of top 30 \npredictions for each model")
  ggsave("timeseries_plots/published_methods/definitive_endoderm_discordance_top_30.pdf", width = 5, height = 7)
  
  
  # perturbation scores versus lit review
  heatmap_data %>% 
    write.csv("timeseries_plots/published_methods/definitive_endoderm_ps_vs_screen_top_30.csv")
  heatmap_data %>% 
    subset(T, select = c("cell_type", "effect_direction")) %>%
    table() %>%
    write.csv("timeseries_plots/published_methods/definitive_endoderm_ps_vs_screen_top_30_unique.csv")
  heatmap_data %>% 
    mutate(prediction_for = paste0("Prediction about:\n", cell_type)) %>%
    subset(predicted_effect!="none") %>%
    subset(predicted_effect!="Not in top 30") %>%
    subset(!is.na(predicted_effect)) %>%
    mutate(accuracy = "Incorrect\n(no effect)") %>%
    mutate(accuracy = ifelse(effect_direction != "none",           "Incorrect\n(wrong direction)", accuracy)) %>%
    mutate(accuracy = ifelse(effect_direction == predicted_effect, "Correct",                     accuracy)) %>%
    mutate(accuracy = forcats::fct_relevel(accuracy, rev(sort(unique(accuracy)))[c(1,3,2)])) %>%
    mutate(method = factor(method, levels = sort(as.character(unique(method))))) %>%
    ggplot() + 
    geom_bar(aes(x=as.character(prediction_timescale), fill = accuracy)) + 
    scale_fill_manual(values = c( "Correct"="green", "Incorrect\n(wrong direction)" = "red", "Incorrect\n(no effect)"="gray")) +
    facet_grid(prediction_for~regression_method, scale = "free") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y="Number of genes", x = "Simulation duration") +
    ggtitle("Literature review of top 30 \npredictions for each model")
  ggsave("timeseries_plots/published_methods/definitive_endoderm_lit_review_top_30.pdf", width = 10, height = 4)
}

# Paul et al. 2015 myelopoiesis data
{
  embeddings = load_embeddings("1.5.1_1")
  embeddings$regression_method %<>% gsub("docker____ekernf01/ggrn_docker_backend_", "", .)
  embeddings$regression_method %<>% gsub("dictys", "dictys_dynamics", .)
  embeddings$regression_method %<>% gsub("sckinetics", "sckinetics_dynamics", .)
  
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
    group_by(super_cell_type, perturbation, regression_method) %>%
    summarize(
      perturbation_score = mean(perturbation_score, na.rm = T), 
      lit_review = unique(Annotation_summary), 
      included_co = unique(Included.in.CellOracle.lit.review.)
    ) %>%   
    ungroup() %>%
    tidyr::pivot_wider(names_from = "super_cell_type", values_from = "perturbation_score") %>% 
    group_by(regression_method) %>%
    mutate(GM_abs_rank = rank(-abs(GM)), ME_abs_rank = rank(-abs(ME)) ) 
  
  # Direction of predicted effect
  prettify = function(x) x %>% gsub("-1", "-", .) %>% gsub("1", "+", .) %>% gsub("0", "none", .)
  genes_in_order = overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>%
    subset(abs(GM)>0|abs(ME)>0) %>% 
    group_by(perturbation) %>%
    summarise(rank_by = abs(sum(sign(GM)))+abs(sum(sign(ME)))) %>%
    dplyr::arrange(rank_by) %>%
    extract2("perturbation")
  overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>%
    subset(abs(GM)>0|abs(ME)>0) %>% 
    mutate(perturbation=factor(perturbation, levels = genes_in_order)) %>%
    mutate(GM = prettify(sign(GM))) %>%
    mutate(ME = prettify(sign(ME))) %>%
    tidyr::pivot_longer(cols = c("GM", "ME"), names_to = "cell_type", values_to = "predicted_effect") %>%
    mutate(method = regression_method %>% gsub("_", " ", .)) %>%
    extract(c("method", "cell_type", "perturbation", "predicted_effect")) %>% 
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect="none")) %>%
    ggplot() + 
    geom_tile(aes(x = method, y = perturbation, fill = predicted_effect)) + 
    scale_fill_manual(values = c("-"="blue", "+"="yellow", "none" = "gray"), na.value = "gray") +
    facet_wrap(~cell_type) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave("timeseries_plots/published_methods/paul_discordance_top_30.pdf", width = 4, height = 7)
  
  # Predictions vs lit review
  overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>% 
    subset(abs(GM)>0|abs(ME)>0) %>% 
    write.csv("timeseries_plots/published_methods/paul_ps_vs_screen_top_30.csv")
  overall_scores %>% 
    subset(abs(GM)>0|abs(ME)>0) %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30, select = c("perturbation", "lit_review", "included_co")) %>%
    distinct() %>%
    extract(c("lit_review", "included_co")) %>%
    table() %>%
    write.csv("timeseries_plots/published_methods/paul_ps_vs_screen_top_30_unique.csv")
  overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>% 
    subset(abs(GM)>0|abs(ME)>0) %>% 
    mutate(predicted_effect_on_GM = GM_abs_rank<=30, predicted_effect_on_ME = ME_abs_rank<=30) %>%
    tidyr::pivot_longer( cols = c("predicted_effect_on_GM", "predicted_effect_on_ME"), names_to = "cell_type", values_to = "predicted_effect" ) %>% 
    subset(predicted_effect) %>%
    mutate(cell_type = gsub("predicted_effect_on_", "", cell_type)) %>%
    mutate(cell_type = paste0("Prediction about:\n", cell_type)) %>%
    extract(c("lit_review", "regression_method", "cell_type")) %>%
    rename(known_role = lit_review) %>%
    mutate(known_role = ifelse(is.na(known_role), "Not reviewed", known_role)) %>%
    table() %>%
    as.data.frame() %>%
    mutate(known_role = forcats::fct_relevel(known_role, "Not reviewed", after = 0)) %>%
    mutate(known_role = forcats::fct_relevel(known_role, "No known role", after = 0)) %>%
    ggplot() + 
    geom_bar(aes(x=regression_method,y=Freq, fill = known_role), stat="identity") + 
    scale_fill_manual(values = c(
      "GM"="green",
      "ME"="red",
      "GM & ME, stemness"="yellow",
      "DC" = "blue",
      "No known role"= "black",
      "Not reviewed" = "gray")
    ) +
    labs(fill="Known roles") +
    facet_grid(cell_type~"", scale = "free_y") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Number of genes") +
    ggtitle("Literature review of top 30 \npredictions for each model")
  ggsave("timeseries_plots/published_methods/paul_lit_review_top_30.pdf", width = 4, height = 6)
}

# =========================== Figure 5: cell type-specific networks ===============================
{
  X = collect_experiments(paste0("1.4.1_", 0:4))
  X[X[["pearson_top_20"]]!=0, ]
  X %<>% add_network_cell_type_metadata
  X %<>% make_the_usual_labels_nice
  X %<>% subset(prediction_timescale==10)
  X$pearson %<>% tidyr::replace_na(0)
  X$pearson_top_20 %<>% tidyr::replace_na(0)
  X$pearson_top_100 %<>% tidyr::replace_na(0)
  X$pearson_top_200 %<>% tidyr::replace_na(0)
  X %<>%
    group_by(prediction_timescale, perturbation_dataset, network_datasets, network_source, network_tissue, network_pretty, network_cell_type_matches ) %>%
    summarise(across(METRICS, mean, na.rm = T))
  X %<>% tidyr::pivot_longer(cols = all_of(METRICS), names_to = "metric")
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



