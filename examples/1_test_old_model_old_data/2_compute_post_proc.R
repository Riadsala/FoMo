library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")

options(mc.cores =4, digits = 2)

############################################################################

datasets <- c("kristjansson2014plos", "tagu2022cog", "hughes2024rsos", "clarke2022qjep") 
models <- c("1_0", "1_1", "1_2")

############################################################################

# wrapper function for computing train/test accuracy
compare_FoMo_accuracy <- function(dataset) {
  
  d <- import_data(dataset)
  
  dout <- tibble()
  
  for (modelver in models)
  {
    
    m <- readRDS(paste0("scratch/", dataset, "_train_", modelver, ".model"))
    t <- readRDS(paste0("scratch/", dataset, "_test_", modelver, ".model"))
    
    pred <- summarise_postpred(list(training = m, testing = t), d, 
                                 get_sim = FALSE, draw_sample_frac = 1)
    
    acc <- compute_acc(pred$acc, compute_hpdi = FALSE) %>% 
      mutate(model = paste("FoMo",  modelver))
    
    rm(m, t, pred) 
    
    dout <- bind_rows(dout, bind_rows(acc, acc) %>%
                        mutate(data = dataset))

  
  }
  
  return(dout)
  
 }

############################################################################

for (ds in datasets) {
  
  d_acc_c2022 <- compare_FoMo_accuracy(ds)
  write_csv(d_acc_c2022, paste0("scratch/post_acc_", ds, ".csv"))
}

############################################################################
# compute simulated run statistics
############################################################################

for (model_ver in models) {
  
  rl <- tibble()
  
  for (ds in datasets) {
    
    # load dataset and model
    d <- import_data(ds)
    m <- readRDS(paste0("scratch/", ds, "_train_", model_ver, ".model"))
    
    # get simulation data for model
    pred <- summarise_postpred(list(training = m, testing = m), d, 
                               get_sim = TRUE, draw_sample_frac = 0.01)
    
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

for (model_ver in models) {
  
  iisv <- tibble()
  
  for (ds in datasets) {
    
    # load dataset and model
    d <- import_data(ds)
    m <- readRDS(paste0("scratch/", ds, "_train_", model_ver, ".model"))
    
    # get simulation data for model
    pred <- summarise_postpred(list(training = m, testing = m), d, 
                               get_sim = TRUE, draw_sample_frac = 0.01) 
    
    # compute empirical run statistics
    iisve <- get_iisv_over_trials(d$found) 
    
    # compute simulated run statistics
    iisvp <- get_iisv_over_trials(pred$sim %>%
                                    # it would be great to remove this line
                                    filter(is.finite(x)))
    
    # bind everything together
    bind_rows(iisve %>%  mutate(x = "human"),
              iisvp  %>% mutate(x = "model")) %>%
      mutate(dataset = ds,
             model_ver = model_ver) %>% 
      bind_rows(iisv) -> iisv
    
  }
  
  write_csv(iisv, paste0("scratch/iisv_statistics", model_ver, ".csv"))
  
}
