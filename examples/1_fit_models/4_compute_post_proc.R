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

# datasets <- c("kristjansson2014plos", "tagu2022cog", "hughes2024rsos", "clarke2022qjep") 
datasets <- "kristjansson2014plos"

############################################################################


extract_and_save_predictions <- function(dataset) {
  
  # wrapper function for computing train/test accuracy for each version
  # of FoMo for a given dataset
  
 
  
  # get list of model versions to compute over
  folder <- paste0("scratch/models/", dataset, "/sim/")
  models <- dir(folder, ".csv")
  models <- get_models_in_dir(folder, mode)
  m <- read_cmdstan_csv(paste0(folder, models))
  
  # tidy
  as_tibble(m$generated_quantities) %>%
    pivot_longer(everything(), names_to = "param") %>%
    separate(param, into = c("draw", "param"), sep = "\\.") %>%
    separate(param, into = c("param", "row"), sep = "\\[") %>%
    mutate(draw = parse_integer(draw),
           row = parse_number(row)) -> genquant
  
  genquant %>% 
    filter(param != "Q") %>%
    pivot_wider(names_from = "param")
  
  
  # create output folder
  outfolder <- paste0("scratch/post/", dataset)
  
  if(!dir.exists("scratch/post/")) {
    dir.create("scratch/post/")
  }
  
  # create save folder if it doesn't yet exist
  if(!dir.exists(outfolder)) {
    dir.create(outfolder)
  }
  
  # read in models and extract post predictions
  for (modelver in models)
  {
    print(paste("... model version ", modelver))
    
 
    dataset <- "kristjansson2014plos"
    
    # get all model predictions
    pred <- extract_pred(dataset, modelver)
    
    pred$dataset <- dataset
    pred$model_ver <- modelver
    
    # save
    print("saving data")
    saveRDS(pred, paste0(outfolder, "/pred_", mode, modelver, ".rds"))
    
    # summarise accuracy and save
    print("summarising accuracy....")
    acc <- summarise_acc(pred)
    write_csv(acc, paste0(outfolder, "/acc_", mode, modelver, ".csv"))
    rm(m, t, pred) 
    
  }
}

compute_iisv_and_run_statistics <- function(dataset){
  
  #####################################################
  # Process empirical data first
  
  d <- import_data(dataset)
  
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
  
  # get list of model versions to compute over
  folder <- paste0("scratch/models/", dataset, "/")
  mode <- "train"
  models <- get_models_in_dir(folder, mode)
  
  # update folder to point to scratch/post
  folder <- paste0("scratch/post/", dataset, "/")
  
  for (modelver in models) {
    
    # get simulation data for model
    pred <- readRDS(paste0(folder, "pred_", mode, modelver, ".rds")) 
    pred <- pred$trialwise %>% filter(.draw == 1)
    
    # compute simulated run statistics
    rlp <- get_run_info_over_trials(pred) %>%
      group_by(person, condition) %>%
      summarise(max_run_length = mean(max_run_length),
                num_runs = mean(n_runs),
                mean_bestr = mean(best_r),
                mean_pao = mean(pao),
                .groups = "drop")
    
    # compute simulated iisv statistics
    #iisvp <- get_iisv_over_trials(pred)
    
    # bind everything together
    rl %>% bind_rows(rlp %>% mutate(z = paste0("v",  modelver))) -> rl
    
  }
  
  rl %>% 
    pivot_longer(c(max_run_length, num_runs, mean_bestr, mean_pao), names_to = "statistic") %>%
    pivot_wider(names_from = z) -> rl
  
  write_csv(rl, paste0(folder, "run_statistics.csv"))
  
  }

############################################################################
# extract model predictions
############################################################################
for (ds in datasets) {
  
  print(paste("Obtaining posterior predictions for dataset ", ds))
  
  # first, extract and save accuracy
  print("***** Computing accuracy *****")
  extract_and_save_predictions(ds)
  
  print("***** Computing iisv and run statistics *****")
  compute_iisv_and_run_statistics(ds)

}

