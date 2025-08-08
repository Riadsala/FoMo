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
datasets <- c( "hughes2024rsos", "clarke2022qjep", "kristjansson2014plos", "tagu2022cog") 
############################################################################

extract_and_save_predictions <- function(dataset, draws_for_sim = 1) {
  
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

  # read in models and extract post predictions
  for (file in models)
  {
    
    # get model version from file name
    modelver <- str_extract(file, "\\d_\\d")
    
    # get kappa from filename
    kappa <- str_extract(file, "_k[\\d]*")
    modelver <- if_else(is.na(kappa), modelver, paste0(modelver, kappa))
    
    print(paste("... model version ", modelver))

    # get all model predictions
    pred <- extract_pred(dataset, modelver, folder)
    
    # save
    print("saving predictions")
    saveRDS(pred, paste0(outfolder, "/pred_", modelver, ".rds"))
  
    rm(pred) 
    
  }
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

