
fit_model <- function(dataset, fomo_ver, mode = "all",
                      model_components = c("spatial", "item_class"),
                      iter = 1000) {
  
  if (class(dataset) != "list") {
    
    d <- import_data(dataset)
    
  } else {
    d <- dataset
    dataset <- "unknown"
  }
  
  # if mode = "all", fit data to everything
  # if mode = "traintest, fit to training set then eval on test set
  
  fomo_ver_str <- str_replace(fomo_ver, "\\.", "_" )
  mod <- cmdstan_model(paste0("../../models/multi_level/FoMo", fomo_ver_str, ".stan"))
  
  # check if we are carrying out a prior model only
  if (fomo_ver_str == "0_0") {
    fxdp = TRUE
  } else {
    fxdp = FALSE
  }
  
  if (mode == "all") {
    
    # check if we have already computed and saved this
    
    d_list <- prep_data_for_stan(d$found, d$stim, model_components)
    d_list <- add_priors_to_d_list(d_list, modelver = fomo_ver)
    
    m <- mod$sample(data = d_list, 
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 10, 
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3,
                    fixed_param = fxdp)
    
    filename_all <- paste0("scratch/", dataset, "_all_", fomo_ver_str, ".model")
    m$save_object(filename_all)
    
  } else if (mode == "traintest") {
    
    d_list <- prep_train_test_data_for_stan(d)
    
    d_list$training <- add_priors_to_d_list(d_list$training, modelver = fomo_ver)
    d_list$testing  <- add_priors_to_d_list(d_list$testing,  modelver = fomo_ver)
    
    # run model
    m_train <- mod$sample(data = d_list$training, 
                          chains = 4, parallel_chains = 4, threads = 4,
                          refresh = 0, 
                          iter_warmup = iter, iter_sampling = iter,
                          sig_figs = 3,
                          fixed_param = fxdp)
    
    m_test <- mod$generate_quantities(m_train, data = d_list$testing, seed = 123)
    
    # save
    filename_train <- paste0("scratch/", dataset, "_train_", fomo_ver_str, ".model")
    m_train$save_object(filename_train)
    
    filename_test <- paste0("scratch/", dataset, "_test_", fomo_ver_str, ".model")
    m_test$save_object(filename_test)
    
  }
  
}
