library(cmdstanr)


fit_model <- function(dataset, fomo_ver, mode = "all",
                      iter = 1000, iter_genquant = 10) {
  
  #######################################################################
  # wrapper function for loading d_list, fitting model, train, test etc. 
  
  # - dataset should be either a label, or a list containing dataset$name
  # - mode = "all" means that training and testing will be be done on all 
  # the data. "split" or "traintest" means that training and testing 
  # take place on separate partitions of the data
  # - iter = number of iterations during model fitting
  # - iter_genquant = number of posterior prediction/simulation samples
  # - required
  #######################################################################
  
  #######################################################################
  # getting setup etc.
  
  # check that if dataset is a list, did we provide a name?
  if (class(dataset) == "list") {
    if (is.null(dataset$name)) {
      print("ERROR: please provide a dataset name in d")
      return()
    }
  }
  
  # get dataset name
  dataset_name <- get_dataset_name(dataset)
  
  # get filepath for the stan files for fitting multi-level model 
  # for for simulating new data
  paths <- get_paths(dataset_name)
  
  # Change 1.x to 1_x if required
  fomo_ver_str <- str_replace(fomo_ver, "\\.", "_" )
  
  # load the two stan files
  mod_fit <- cmdstan_model(paste0(paths$model, "FoMo", fomo_ver_str, ".stan"))
  mod_sim <- cmdstan_model(paste0(paths$simul, "FoMo", fomo_ver_str, ".stan"))
  
  # check if we are carrying out a prior model only
  if (fomo_ver_str == "0_0") {
    fxdp = TRUE
  } else {
    fxdp = FALSE
  }
  
  #######################################################################
  # load data and fit model
  
  # load the pre-computed d_list and add required priors
  d_list <- get_list(dataset, mode, "training")
  d_list <- add_priors_to_d_list(d_list, modelver = fomo_ver, model_path = paths$model)
  
  # fit the model
  m <- mod_fit$sample(data = d_list,
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10,
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3,
                  fixed_param = fxdp,
                  output_dir = paths$out_fit,
                  output_basename = paste(dataset_name, fomo_ver_str, sep=""))
  
  # now save
  m$save_object(paste0(paths$out_fit, fomo_ver_str, ".model"))
  
  ###########################################################################
  # now create generated quantities from fitted model

  # do we need to load in a new d_list for testing?
  if (mode %in%  c("split", "traintest")) {
    d_list <- get_list(dataset, mode, "testing")
  }
  
  # randomly sample some draws to calculate generated quantities for
  draws_matrix <- posterior::as_draws_matrix(m$draws())
  idx <- sample(nrow(draws_matrix), iter_genquant)
    
  p <- mod_sim$generate_quantities(fitted_params = draws_matrix[idx,], 
                                        data = d_list, 
                                        seed = 123,
                                        output_dir = paths$out_sim,
                                        output_basename = paste(dataset_name, fomo_ver_str, sep=""))
    
  p$save_object(paste0(paths$out_sim, dataset_name, fomo_ver_str, ".model"))

}

get_paths <- function(ds) {
  
  # get model directory path
  model_path <- "../"
  simul_path <- "../"
  
  for (ii in 1:3) {
    
    if ("models" %in% dir(model_path)) {
      model_path <- paste0(model_path, "models/")
      
      break
    } else {
      model_path <- paste0(model_path, "../")
    }
  }
  
  simul_path <- paste0(model_path, "simulate/")
  model_path <- paste0(model_path, "multi_level/")
  
  # now sort out output path
  outpt_path <- paste0("scratch/models/", ds)
  # create save folder if it doesn't yet exist

    dir.create(outpt_path)
    dir.create(paste0(outpt_path, "/fit/"))
    dir.create(paste0(outpt_path, "/sim/"))

  
  return(list(model = model_path, 
              simul = simul_path,
              out_fit = paste0(outpt_path, "/fit/"),
              out_sim = paste0(outpt_path, "/sim/")))
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
