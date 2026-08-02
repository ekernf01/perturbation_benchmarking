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

# =========================== Figure 2, perturbation scores and lit reviews for published methods ==========
dir.create("timeseries_plots/published_methods")
# Both datasets; quantitative evals
{
  X = collect_experiments(paste0("1.5.1_", 0:2)) %>% make_the_usual_labels_nice
  X %<>%
    group_by(regression_method, perturbation_dataset, prediction_timescale) %>%
    summarise(across(METRICS, mean, na.rm = T))
  X %<>% tidyr::pivot_longer(cols = all_of(METRICS), names_to = "metric")
  X[["metric"]] %<>% gsub("_", " ", .)
  X[["metric"]] %<>% factor(levels = gsub("_", " ", METRICS))
  plot = ggplot(X) +
    geom_bar(aes(x = gsub("_", " ", regression_method), y = value, fill = factor(as.character(prediction_timescale), levels = c(0:5, 10)), group = prediction_timescale), position = "dodge", stat = "identity") +
    facet_grid(~perturbation_dataset, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))  +
    ylab("Pearson correlation between\npredicted and observed log FC") + 
    xlab("") +
    geom_hline(yintercept=0) + 
    labs(fill="Simulation\nduration") + 
    scale_fill_manual(values = viridis::viridis(6), guide = "legend")
  print(plot)
  ggsave('timeseries_plots/published_methods/pearson.svg', plot, width = 6, height = 3)
}

# Definitive Endoderm, effects on cell type 
{
  embeddings = load_embeddings("1.5.1_0") %>% make_the_usual_labels_nice()
  relevant_columns = c("cell_type",
                       "perturbation", 
                       "prediction_timescale",
                       "regression_method",
                       "perturbation_score",
                       "brunello.neg.lfc",
                       "brunello.pos.lfc",
                       "gecko.neg.lfc",
                       "gecko.pos.lfc"
  )
  embeddings = as.data.frame(embeddings)[relevant_columns]
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
      Included.in.literature.review. = unique(Included.in.literature.review), 
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
  ggsave("timeseries_plots/published_methods/definitive_endoderm_ps_vs_screen.pdf", width = 8, height = 3)
  
  # perturbation scores overall distribution
  overall_scores %>% 
    ggplot() + 
    stat_ecdf(aes(x = perturbation_score, color = factor(as.character(prediction_timescale), levels = c(1:5, 10)), group = prediction_timescale)) + 
    facet_grid(gsub("_| ", "\n", regression_method)~cell_type) + 
    labs(y = "Cumulative density") + 
    geom_vline(xintercept = 0, color = "red") + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
    labs(color="Simulation\nduration") + 
    scale_color_manual(values = viridis::viridis(6), guide = "legend")
  ggsave("timeseries_plots/published_methods/definitive_endoderm_ps_distribution.pdf", width = 4, height = 8)
  
  
  # Compare top perturbation scores versus each other
  # First, add an extra baseline from differential expression.
  # This fit awkwardly into the code because there is no gene expression prediction involved, so the normal eval code doesn't apply.
  {
    differential_expression = rbind(
      read.csv("timeseries_plots/top30_differential_expression_endoderm.csv")    %>% mutate(cell_type = "endoderm", perturbation_score_rank = rank(logfoldchanges)),
      read.csv("timeseries_plots/top30_differential_expression_mesendoderm.csv") %>% mutate(cell_type = "mesendoderm", perturbation_score_rank = rank(logfoldchanges))
    ) %>% rename( 
      perturbation =  names, 
    ) %>% mutate(
      perturbation_score = -logfoldchanges,
      prediction_timescale = 1, 
      regression_method = "differential expression", 
      top_30 = "yes",
      li_et_al_biggest_effect = NA
    ) %>% dplyr::select(
      cell_type, perturbation, perturbation_score, prediction_timescale, regression_method, li_et_al_biggest_effect, top_30, perturbation_score_rank
    )
    lit_review = read.csv("../../perturbation_data/perturbations/definitive_endoderm/lit_review.csv")
    lit_review %<>% rename(cell_type = Cell.type.affected)
    lit_review %<>% mutate(cell_type = gsub("primitive streak", "mesendoderm", tolower(cell_type)))
    lit_review %<>% tidyr::complete(perturbation, cell_type, fill = list(effect_direction="none"))
    differential_expression %<>% merge(lit_review, by = c("perturbation", "cell_type"), all.x = T)
    differential_expression %<>% mutate(Included.in.literature.review. = ifelse(is.na(Included.in.literature.review), "no", Included.in.literature.review))
    overall_scores = rbind(overall_scores, differential_expression[colnames(overall_scores)])
  }
  
  # Now, some code to make the heatmap and dress it up. 
  prettify = function(x) x %>% gsub("-1", "LoF inhibits differentiation", .) %>% gsub("1", "LoF promotes differentiation", .) %>% gsub("0", "none", .)
  # Ordering heatmaps is always a nightmare.
  # This makes a rough ordering according to th lit review results. 
  heatmap_data = overall_scores %>% 
    subset(perturbation_score_rank<=30) %>% 
    mutate(predicted_effect = prettify(sign(perturbation_score))) %>%
    mutate(method = regression_method %>% gsub("_", " ", .)) %>%
    mutate(method = interaction(regression_method, prediction_timescale)) %>%
    extract(c("method", "regression_method", "prediction_timescale", "cell_type", "perturbation", "predicted_effect")) %>% 
    ungroup() %>%
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect="none")) %>%
    merge(lit_review, by = c("cell_type", "perturbation"), all.x = T, all.y = F)
  genes_in_order = heatmap_data %>% 
    dplyr::arrange(effect_direction) %>%
    extract2("perturbation") %>%
    unique()
  
  simplified_prediction = c(
    "none" = 0,
    "Not in top 30" = 0, 
    "GoF promotes differentiation" = -1,
    "LoF inhibits differentiation" = -1,
    "LoF affects expression state" = -1,
    "LoF promotes differentiation" = 1
  )
  heatmap_data %<>% 
    mutate(perturbation=factor(perturbation, levels = genes_in_order)) %>%
    mutate(predicted_effect = gsub("none", "Not in top 30", predicted_effect)) %>%
    mutate(effect_simple = simplified_prediction[predicted_effect])
  
  simplified_plot_data = 
    heatmap_data[c("method", "effect_simple", "cell_type", "perturbation")] %>%
    tidyr::pivot_wider(id_cols = c("cell_type", "perturbation"), names_from = "method", values_from = "effect_simple") %>%
    extract(-(1:2)) %>%
    cor %>%
    as.data.frame() %>%
    tibble::rownames_to_column() %>%
    tidyr::pivot_longer(cols = -rowname) %>%
    rename(method1 = rowname, method2 = name, correlation = value) 
  ggplot(simplified_plot_data) + 
    geom_tile(aes(x = method1, y = method2, fill = correlation)) + 
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0
    ) + 
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) + 
    ggtitle("Correlation in prediction \ndirection across methods") + 
    coord_fixed()
    
  ggsave("timeseries_plots/published_methods/definitive_endoderm_discordance_top_30.pdf", width = 6, height = 6)
  
  
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

# Paul et al. 2015 myelopoiesis data; effect on cell type
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
  
  # Compare top perturbation scores versus each other
  # First, add an extra baseline from differential expression.
  # This fit awkwardly into the code because there is no gene expression prediction involved, so the normal eval code doesn't apply.
  {
    differential_expression = rbind(
      read.csv("timeseries_plots/top30_differential_expression_paul_ME.csv") %>% mutate(super_cell_type = "ME"),
      read.csv("timeseries_plots/top30_differential_expression_paul_GM.csv") %>% mutate(super_cell_type = "GM"),
      read.csv("timeseries_plots/top30_differential_expression_paul_DC.csv") %>% mutate(super_cell_type = "DC")
    ) %>% rename( 
      perturbation =  names, 
    ) %>% mutate(
      perturbation_score = -logfoldchanges,
      regression_method = "differential expression", 
    ) %>% dplyr::select(
      super_cell_type, perturbation, perturbation_score, regression_method
    ) %>% tidyr::pivot_wider(
      values_from = "perturbation_score",
      names_from = "super_cell_type", 
      id_cols = c("perturbation", "regression_method")
    )
    lit_review = read.csv("../../perturbation_data/perturbations/paul1/screen.csv")
    differential_expression %<>% merge(lit_review, by = "perturbation", all.x = T, all.y = F)
    differential_expression %<>% mutate(lit_review = Annotation_summary, included_co = Included.in.CellOracle.lit.review.)
    differential_expression %<>% select(c(
      "perturbation" ,     "regression_method" , "lit_review"   ,     "included_co" ,      "GM"      ,          "ME"                   
      ))
    # This gets added into "overall_scores" shortly below. 
  }
  
  
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
    rbind(differential_expression) %>% # This is where we add in the differential expression baseline
    mutate(GM_abs_rank = rank(-abs(GM)), ME_abs_rank = rank(-abs(ME)) ) 
  
  # Direction of predicted effect
  genes_in_order = overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>%
    subset(abs(GM)>0|abs(ME)>0) %>% 
    group_by(perturbation) %>%
    summarise(rank_by = abs(sum(sign(GM)))+abs(sum(sign(ME)))) %>%
    dplyr::arrange(rank_by) %>%
    extract2("perturbation")
  # TODO: SIMPLIFY THIS GRAPHIC
  correlations = overall_scores %>% 
    subset(GM_abs_rank<=30|ME_abs_rank<=30) %>%
    subset(abs(GM)>0|abs(ME)>0) %>% 
    mutate(perturbation=factor(perturbation, levels = genes_in_order)) %>%
    mutate(GM = (sign(GM))) %>%
    mutate(ME = (sign(ME))) %>%
    tidyr::pivot_longer(cols = c("GM", "ME"), names_to = "cell_type", values_to = "predicted_effect") %>%
    mutate(method = regression_method %>% gsub("_", " ", .)) %>%
    extract(c("method", "cell_type", "perturbation", "predicted_effect")) %>% 
    tidyr::complete(method, perturbation, cell_type, fill = list(predicted_effect=0)) %>%
    unique() %>%
    tidyr::pivot_wider(id_cols = c("cell_type","perturbation"), values_from = "predicted_effect", names_from = "method") 
  correlations = rbind(
    cbind(tibble::rownames_to_column(data.frame(cor(subset(correlations, cell_type == "GM")[-(1:2)]))), cell_type = "GM"),
    cbind(tibble::rownames_to_column(data.frame(cor(subset(correlations, cell_type == "ME")[-(1:2)]))), cell_type = "ME")
  ) %>% 
    tidyr::pivot_longer(-c(cell_type, rowname), values_to = "correlation") %>%
    mutate(name=gsub("\\.", " ", name)) %>% 
    subset(name!=rowname)
  correlations %>%
    ggplot() + 
    geom_tile(aes(x = rowname, y = name, fill = correlation)) + 
    facet_wrap(~cell_type) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
    xlab("") + ylab("")  + 
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0
    ) + 
    coord_fixed() +
    ggtitle("Correlation in prediction direction\nacross methods (Paul data)")
  ggsave("timeseries_plots/published_methods/paul_discordance_top_30.pdf", width = 7, height = 6)
  
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
    mutate(known_role = forcats::fct_relevel(known_role, "no known effect", after = 0)) %>%
    ggplot() + 
    geom_bar(aes(x=regression_method %>% gsub("_", " ", .), y=Freq, fill = known_role), stat="identity") + 
    scale_fill_manual(values = c(
      "GM & ME"       = "yellow",
      "GM & DC"       = "cyan",
      "GM & ME & DC" = "brown",
      "GM"="green",
      "ME"="red",
      "Stemness" = "orange",
      "DC" = "blue",
      "ME & DC" = "purple",
      "no known effect"= "black",
      "unclear"= "gray",
      "Not reviewed" = "gray")
    ) +
    labs(fill="Known roles") +
    facet_grid(cell_type~"", scale = "free_y") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ylab("Number of genes") +
    xlab("") +
    ggtitle("Literature review of top 30 \npredictions for each model")
  ggsave("timeseries_plots/published_methods/paul_lit_review_top_30.pdf", width = 6, height = 6)
}

