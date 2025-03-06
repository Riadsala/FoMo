library(tidyverse)
library(cmdstanr)
library(loo)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/fit_model.R")

# read in data
datasets <- "ac_test"

extract_and_save_predictions <- function(dataset) {
  # wrapper function for computing train/test accuracy for each version
  # of FoMo for a given dataset
  
  d <- readRDS(paste0("scratch/data/", dataset, ".RDS"))
  
  # get list of model versions to compute over
  folder <- paste0("scratch/models/", dataset, "/")
  mode <- "train"
  models <- unlist(dir(folder))
  models <- models[str_detect(models, mode)]
  models <- str_extract(models, "1_[0-9]")
  
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
    
    m <- readRDS(paste0("scratch/models/", dataset, "/train", modelver, ".model"))
    t <- readRDS(paste0("scratch/models/", dataset, "/test", modelver, ".model"))
    
    # get all model predictions
    pred <- extract_pred(list(training = m, testing = t), d)
    
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


############################################################################
# extract model predictions
############################################################################
for (ds in datasets) {
  
  print(paste("Obtaining posterior predictions for dataset ", ds))
  extract_and_save_predictions(ds)
  
}













