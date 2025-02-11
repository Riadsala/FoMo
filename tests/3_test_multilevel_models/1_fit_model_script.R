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

# create scratch folder if it does not yet exist
if (!file.exists("scratch"))  dir.create("scratch")


# create parameters for simulation
exp <- tibble(
  n_people = 8,
  n_conditions = 1,
  condition_name = "simple test",
  n_trials_per_cond = 2,
  n_item_class = 2, 
  n_item_per_class = 5
)

foraging <- tribble(
  ~param, ~mu, ~sd,
  "bA", c(0.7, 0.3, 0, 0), 0.2,
  "bS", 1, 1,
  "rho_delta", 20, 2,
  "rho_psi", 1, 1
)

abs_dir_tuning = list(kappa = rep(10, 4), theta = rep(1, 4))

# initial bias params
inital_sel_params <- tibble(
  a1x = 2,
  b1x = 2,
  a2x = 1,
  b2x = 10,
  a1y = 2,
  b1y = 2,
  a2y = 10,
  b2y = 1) 

params <- list(exp = exp, 
               foragaging = foraging, 
               abs_dir_tuning = abs_dir_tuning, 
               inital_sel_params = inital_sel_params)

# now simulate data
d <- sim_foraging_people(params,
                         rel_proximity = FALSE,
                         filename = "test_anna") 

# fit all the models
m <- fit_model(d, fomo_ver = "1.0", mode = "all", iter = 500) 
m <- fit_model(d, fomo_ver = "1.0", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.1", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.2", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.3", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.4", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.5", mode = "traintest", iter = 500) 
