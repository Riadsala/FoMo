
fit_model <- function(dataset, fomo_ver, mode = "all",
                      model_components = c("spatial", "item_class"),
                      iter = 1000) {
  
  # check that if dataset is a list, did we provide a name?
  if (class(dataset) == "list") {
    if (is.null(dataset$name)) {
      print("ERROR: please provide a dataset name in d")
      return()
    }
  }
  
  # get dataset name
  dataset_name <- get_dataset_name(dataset)
  
  # create save folder if it doesn't yet exist
  if(!dir.exists(paste0("scratch/models/", dataset_name))) {
    dir.create(paste0("scratch/models/", dataset_name))
  }
  # get model directory path
  model_path <- "../"
  
  for (ii in 1:3) {
    
    if ("models" %in% dir(model_path)) {
      model_path <- paste0(model_path, "models/")
      break
    } else {
      model_path <- paste0(model_path, "../")
    }
  }
  
  # load in the Stan model
  fomo_ver_str <- str_replace(fomo_ver, "\\.", "_" )
  mod <- cmdstan_model(paste0(model_path, "multi_level/FoMo", fomo_ver_str, ".stan"))
  
  # check if we are carrying out a prior model only
  if (fomo_ver_str == "0_0") {
    fxdp = TRUE
  } else {
    fxdp = FALSE
  }
  
  ###########################################################################
  # fit model to data 
  # either "all" or "training"
  
  # get d_list
  d_list <- get_list(dataset, mode, "training")
 
  # add priors to d_list
  d_list <- add_priors_to_d_list(d_list, modelver = fomo_ver, model_path = model_path)

  m <- mod$sample(data = d_list, 
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 10, 
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3,
                    fixed_param = fxdp)
  
  # now save
  if (mode == "all") {
    filename <- paste0("scratch/models/", dataset_name, "/all", fomo_ver_str, ".model")
  } else {
    filename <- paste0("scratch/models/", dataset_name, "/train", fomo_ver_str, ".model")
  }
  
  m$save_object(filename)
  
  ###########################################################################
  # If we are in train-test model, we should now 
  # compute generated quantities for the test data
  ###########################################################################
  if (mode == "traintest") {
    
    d_list <- get_list(dataset, mode, "testing")
    # although we aren't using the priors, the model still
    # expects them to be in the input
    d_list  <- add_priors_to_d_list(d_list, modelver = fomo_ver, model_path = model_path)
    
    m_test <- mod$generate_quantities(m, data = d_list, seed = 123)
    
    filename <- paste0("scratch/models/", dataset_name, "/test", fomo_ver_str, ".model")
    m_test$save_object(filename)

  }
  
  # return model
  if (mode == "all") {
    return(m)
  } else {
    return(list(training = m, testing = m_test))
  }
}

get_list <- function(dataset, mode, stage) {
  
  if (class(dataset) == "character") {
    
    # if we provide a dataset label, load a precomputed d_list
    dlist_folder <- paste0("scratch/d_list/", dataset, "/")
    
    if (mode == "all") {
      d_list <- readRDS(paste0(dlist_folder, "all.rds"))
    } else {
      
      # load either training or testing dlist
      if (stage == "training") {
        d_list <- readRDS(paste0(dlist_folder, "train.rds"))
      } else {
        d_list <- readRDS(paste0(dlist_folder, "test.rds"))
      }
    }
    
  } else {
    # otherwise... compute d_list from scratch
    
    if (mode == "all") {
      d_list <- prep_data_for_stan(d)
    } else {
      d_list <- prep_train_test_data_for_stan(d)
      
      if (stage == "training") {
        d_list = d_list$training
      } else {
        d_list = d_list$testing
      }
    }
  }
}

get_dataset_name <- function(ds) {
  
  # as ds can either be a string (pointing to a dataset name) or
  # an actual list-dataset, we better sort that out
  
  if (class(ds) == "list") {
    ds_name <- ds$name
    
  } else {
    ds_name <- ds 
  }
  
  return(ds_name)
  
}
