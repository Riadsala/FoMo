library(tidyverse)
library(cmdstanr)

# This script reads evaluates model accuracy and computes summaries
# allowing us to compare human and model run statistics and inter-
# item selection vectors. 

source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/post_functions.R")

options(mc.cores = 4, digits = 2)

draws_for_sim <- 2

############################################################################

datasets <- c("hughes2024rsos", "tagu2022cog", "kristjansson2014plos")  #"clarke2022qjep", "hughes2024rsos", "tagu2022cog",

############################################################################

compute_summary_stats <- function(dataset, draws_for_sim = 1) {
  
  # wrapper function for computing train/test accuracy for each version
  # of FoMo for a given dataset
  
  ############################################################################
  # get set up
  
  # get list of model versions to compute over
  folder <- paste0("scratch/models/", dataset, "/sim/")
  models <- dir(folder, ".csv")
  
  # create output folder
  outfolder <- paste0("scratch/post/", dataset)
  
  if(!dir.exists("scratch/post/")) {
    dir.create("scratch/post/")
  }
  
  # create save folder if it doesn't yet exist
  if(!dir.exists(outfolder)) {
    dir.create(outfolder)
  }
  
  ############################################################################
  # compute statistics for empirical data first
  
  d <- import_data(dataset)
  
  # we only want to calcualte these on the test data
  d <- get_train_test_split(d)
  d <- d$testing
  
  # compute empirical run statistics
  rl <- get_run_info_over_trials(d$found) %>%
    group_by(person, condition) %>%
    summarise(max_run_length = mean(max_run_length),
              num_runs = mean(n_runs),
              mean_bestr = mean(best_r),
              mean_pao = mean(pao),
              .groups = "drop") %>%
    pivot_longer(-c(person, condition), names_to = "statistic", values_to = "observed") 

  # compute empirical run statistics
  iisv <- get_iisv_over_trials(d$found) %>%
    mutate(z = "observed")

  
  # tidy up
  rm(d)
  
  # read in models and extract post predictions

  for (file in models)
  {
    
    # get model version from file name
    modelver <- str_extract(file, "\\d_\\d")
    
    # get kappa from filename
    kappa <- str_extract(file, "_k[\\d]*")
    modelver <- if_else(is.na(kappa), modelver, paste0(modelver, kappa))
    
    # load
    print(paste("loading predictions for model version ", modelver))
    pred <- readRDS( paste0(outfolder, "/pred_", modelver, ".rds"))
    
    # compute simulated run statistics
    print("computing predicted run statistics")
    rlp <- get_run_info_over_trials(pred$trialwise %>% filter(.draw < (draws_for_sim+1))) %>%
      group_by(.draw, person, condition) %>%
      summarise(max_run_length = mean(max_run_length),
                num_runs = mean(n_runs),
                mean_bestr = mean(best_r),
                mean_pao = mean(pao),
                .groups = "drop") %>% 
      pivot_longer(-c(.draw, person, condition), 
                   names_to = "statistic", 
                   values_to =  paste0("v",  modelver))
    
    rl <- full_join(rl, rlp)
    
    print("repeat, for fixed first selected...")
    rlp <- get_run_info_over_trials(pred$trialwise_firstfixed %>% filter(.draw < (draws_for_sim+1))) %>%
      group_by(.draw, person, condition) %>%
      summarise(max_run_length = mean(max_run_length),
                num_runs = mean(n_runs),
                mean_bestr = mean(best_r),
                mean_pao = mean(pao),
                .groups = "drop") %>% 
      pivot_longer(-c(.draw, person, condition),
                   names_to = "statistic", 
                   values_to =  paste0("f",  modelver))
    
    rl <- full_join(rl, rlp)
    # 
    # compute empirical run statistics
    iisvp <- get_iisv_over_trials(pred$trialwise %>% filter(.draw <  (draws_for_sim+1))) %>%
      mutate(model_version = paste0("v",  modelver),
             z = "predicted")
    
    iisfp <- get_iisv_over_trials(pred$trialwise_firstfixed %>% filter(.draw <  (draws_for_sim+1))) %>%
      mutate(model_version = paste0("f",  modelver),
             z = "predicted")

    iisv %>%
      bind_rows(iisvp) -> iisv
    
    iisv %>%
      bind_rows(iisfp) -> iisv
    
    rm(pred) 
      
  }
  
  # tidy up run statistics model

  
  # round iisv to 3dp
  iisv %>% mutate(x = round(x, 3), 
                  y = round(y, 3),
                  d2 = round(d2, 3), 
                  theta = round(theta, 3), 
                  psi = round(psi, 3)) -> iisv
  
  # save run statistics
  write_csv(rl,   paste0(outfolder, "/run_statistics.csv"))
  write_csv(iisv, paste0(outfolder, "/iisv_statistics.csv"))
  
}

############################################################################
# extract model predictions
############################################################################
for (ds in datasets) {
  
  print(paste("Obtaining posterior predictions for dataset ", ds))
  compute_summary_stats(ds, draws_for_sim = draws_for_sim)

}

