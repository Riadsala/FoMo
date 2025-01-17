
fit_model <- function(dataset, fomo_ver, mode = "all",
                      model_components = c("spatial", "item_class"),
                      iter = 1000) {
  
  # first, load in the Stan model
  fomo_ver_str <- str_replace(fomo_ver, "\\.", "_" )
  mod <- cmdstan_model(paste0("../../models/multi_level/FoMo", fomo_ver_str, ".stan"))
  
  # check if we are carrying out a prior model only
  if (fomo_ver_str == "0_0") {
    fxdp = TRUE
  } else {
    fxdp = FALSE
  }
  
  ###########################################################################
  # fit model to data 
  # either "all" or "training"
  
  # load d_list
  dlist_folder <- paste0("scratch/d_list/", dataset, "/")
  
  if (mode == "all") {
    d_list <- readRDS(paste0(dlist_folder, "all.rds"))
  } else {
    d_list <- readRDS(paste0(dlist_folder, "train.rds"))
  }

  # add priors to d_list
  d_list <- add_priors_to_d_list(d_list, modelver = fomo_ver)
  

  m <- mod$sample(data = d_list, 
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 10, 
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3,
                    fixed_param = fxdp)
  
  # now save
  dir.create("scratch/models")
  filename <- paste0("scratch/models", dataset, mode, fomo_ver_str, ".model")
  m$save_object(filename)
    
  if (mode == "traintest") {
    
    # now get generated quantities for test data
    d_list <- readRDS(paste0(dlist_folder, "test.rds"))
    #d_list$testing  <- add_priors_to_d_list(d_list$testing,  modelver = fomo_ver)
    
    m_test <- mod$generate_quantities(m_train, data = d_list$testing, seed = 123)
    
    filename <- paste0("scratch/models", dataset, "test", fomo_ver_str, ".model")
    m_test$save_object(filename)

  }
  
}
