library(tidyverse)
library(cmdstanr)

# This script reads evaluates model accuracy and computes summaries
# allowing us to compare human and model run statistics and inter-
# item selection vectors. 

source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/post_functions.R")

options(mc.cores = 4, digits = 2)

############################################################################
datasets <- c("kristjansson2014plos", "tagu2022cog", "hughes2024rsos", "bhat2025",  "clarke2022qjep" ) 
############################################################################

extract_and_save_predictions <- function(dataset) {
  
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
    mutate(z = "observed")
  
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
    print(paste("... model version ", modelver))

    # get all model predictions
    pred <- extract_pred(dataset, modelver, folder)
    
    # save
    print("saving predictions")
    saveRDS(pred, paste0(outfolder, "/pred_", modelver, ".rds"))

    
    # compute simulated run statistics
    print("computing predicted run statistics")
    rlp <- get_run_info_over_trials(pred$trialwise %>% filter(.draw == 1)) %>%
      group_by(.draw, person, condition) %>%
      summarise(max_run_length = mean(max_run_length),
                num_runs = mean(n_runs),
                mean_bestr = mean(best_r),
                mean_pao = mean(pao),
                .groups = "drop") %>% 
      mutate(z = paste0("v",  modelver))
    
    rl %>% bind_rows(rlp) -> rl
    
    print("repeat, for fixed first selected...")
    rlp <- get_run_info_over_trials(pred$trialwise_firstfixed %>% filter(.draw == 1)) %>%
      group_by(.draw, person, condition) %>%
      summarise(max_run_length = mean(max_run_length),
                num_runs = mean(n_runs),
                mean_bestr = mean(best_r),
                mean_pao = mean(pao),
                .groups = "drop") %>% 
      mutate(z = paste0("f",  modelver))
    
    rl %>% bind_rows(rlp) -> rl
    
    # compute empirical run statistics
    iisvp <- get_iisv_over_trials(pred$trialwise %>% filter(.draw == 1)) %>%
      mutate(z = paste0("v",  modelver))
    
    iisv %>% bind_rows(iisvp) -> iisv
    
    rm(pred) 
    
  }
  
  # tidy up run statistics model
  rl %>%  
    select(-.draw) %>%
    pivot_longer(c(max_run_length, num_runs, mean_bestr, mean_pao), names_to = "statistic") %>%
    pivot_wider(names_from = z) -> rl
  
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
  
  # first, extract and save accuracy
  print("***** Extracting preditions *****")
  extract_and_save_predictions(ds)

}

