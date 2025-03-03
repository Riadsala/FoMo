library(tidyverse)


# This script reads evaluates model accuracy and computes summaries
# allowing us to compare human and model run statistics and inter-
# item selection vectors. 

source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/post_functions.R")

options(mc.cores = 4, digits = 2)

############################################################################

# datasets <- c("kristjansson2014plos", "tagu2022cog", "hughes2024rsos", "clarke2022qjep") 
datasets <-"clarke2022qjep"

############################################################################

extract_and_save_predictions <- function(dataset) {
  # wrapper function for computing train/test accuracy for each version
  # of FoMo for a given dataset
  
  d <- import_data(dataset)
  
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
