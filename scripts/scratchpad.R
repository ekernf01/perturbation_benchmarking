library("magrittr")
setwd("~/Desktop/jhu/research/projects/perturbation_prediction/cell_type_knowledge_transfer/")
getRegulators = function(dir){
  read.table(file.path(dir, "genes.regulate.genes")) %>%
    extract2("V4")
}  
cellTypeDirs = list.files("networks/encode_dnase/networks.v12032013/buffer.5000.mm9-120313/", full.names = T)
cellTypeDirs %>% 
  lapply(getRegulators) %>%
  Reduce(union, .) %>%
  unique %>% 
  sort %>%
  is.element(el = c("FOXA2", "SMAD2", "SMAD4", "FOXH1"))
           