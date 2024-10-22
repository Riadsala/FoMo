library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

###################
#


compare_FoMo_accuracy <- function(dataset) {
  
  d <- import_data(dataset)

  m10 <- readRDS(paste0("scratch/", dataset, "_train_1_0.model"))
  t10 <- readRDS(paste0("scratch/", dataset, "_test_1_0.model"))
  
  pred10 <- summarise_postpred(list(training = m10, testing = t10), d, 
                               get_sim = FALSE, draw_sample_frac=0.25)
  acc10 <- compute_acc(pred10$acc) %>% mutate(model = "FoMo 1.0")
  
  rm(m10, t10, pred10) 
  
  m11 <- readRDS(paste0("scratch/", dataset, "_train_1_1.model"))
  t11 <- readRDS(paste0("scratch/", dataset, "_test_1_1.model"))
   
  pred11 <- summarise_postpred(list(training = m11, testing = t11), d, 
                               get_sim = FALSE, draw_sample_frac=0.25)
  acc11 <- compute_acc(pred11$acc) %>% mutate(model = "FoMo 1.1")
  
  rm(m11, t11, pred11)
  
  return(bind_rows(acc10, acc11) %>%
           mutate(data = dataset))
  
 }


d_acc_c2022 <- compare_FoMo_accuracy("clarke2022qjep")