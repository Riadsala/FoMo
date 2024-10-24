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

############################################################################
datasets <- c("kristjansson2014plos", "tagu2022cog", "clarke2022qjep")
############################################################################

# wrapper function for computing train/test accuracy
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

############################################################################

for (ds in datasets)
{
  d_acc_c2022 <- compare_FoMo_accuracy(ds)
  write_csv(d_acc_c2022, paste0("scratch/post_acc_", ds, ".csv"))
}

############################################################################
# compute simulated run statistics
############################################################################

for (model_ver in c("1_0", "1_1")) {
  
  rl <- tibble()
  
  for (ds in datasets) {
    
    # load dataset and model
    d <- import_data(ds)
    m <- readRDS(paste0("scratch/", ds, "_train_", model_ver, ".model"))
    
    # get simulation data for model
    pred <- summarise_postpred(list(training = m, testing = m), d, 
                               get_sim = TRUE, draw_sample_frac = 0.001)
    
    # compute empirical run statistics
    rle <- get_run_info_over_trials(d$found) %>%
      group_by(person, condition) %>%
      summarise(max_run_length = mean(max_run_length))
    
    # compute simulated run statistics
    rlp <- get_run_info_over_trials(pred$sim) %>%
      group_by(person, condition) %>%
      summarise(max_run_length = mean(max_run_length))
    
    # bind everything together
    bind_rows(rle %>% mutate(x = "observed"),
              rlp %>% mutate(x = "predicted")) %>%
      pivot_wider(names_from = "x", values_from = "max_run_length") %>%
      mutate(dataset = ds) %>% 
      bind_rows(rl) -> rl
  }
  
  write_csv(rl, paste0("scratch/run_statistics", model_ver, ".csv"))

}

############################################################################
# compute simulated iisv statistics
############################################################################

for (model_ver in c("1_0", "1_1")) {
  
  iisv <- tibble()
  
  for (ds in datasets) {
    
    # load dataset and model
    d <- import_data(ds)
    m <- readRDS(paste0("scratch/", ds, "_train_", model_ver, ".model"))
    
    # get simulation data for model
    pred <- summarise_postpred(list(training = m, testing = m), d, 
                               get_sim = TRUE, draw_sample_frac = 0.001)
    
    # compute empirical run statistics
    iisve <- get_iisv_over_trials(d$found) %>%
      group_by(condition, found) %>%
      median_hdci(d2)
    
    # compute simulated run statistics
    iisvp <- get_iisv_over_trials(pred$sim) %>%
      group_by(condition, found) %>%
      median_hdci(d2, .width = c(0.53, 0.97))
    
    # bind everything together
    bind_rows(iisve %>% mutate(x = "human"),
              iisvp %>% mutate(x = "model")) %>%
      mutate(dataset = ds) %>% 
      bind_rows(iisv) -> iisv
  }
  
  write_csv(iisv, paste0("scratch/iisv_statistics", model_ver, ".csv"))
  
}
