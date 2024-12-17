# This is a single-use script to help me transfer literature review notes to a better permanent home.

# Endoderm
lit_review = readxl::read_excel("timeseries_figures/definitive_endoderm_ps_vs_screen_top_30_manually_annotated.xlsx")
lit_review = lit_review[c("perturbation", "Included in literature review?", "Cell type affected", "PMID", "Notes")] %>% 
  dplyr::distinct() %>%
  write.csv("../../perturbation_data/perturbations/definitive_endoderm/lit_review.csv")

# TO DO: blood