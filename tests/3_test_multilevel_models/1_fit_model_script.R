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

item_class_weights = list(c(0.7, 0.3, 0, 0))

b_stick = 1

rho_delta = 20 #1 - set to 1 if doing relative proximity = TRUE, 20 if relative proximity = FALSE
sd_rho_delta = 1

rho_psi = 2

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

d <- sim_foraging_people(n_people = 8,
                         n_conditions = 1,
                         cond_lab = c("simple test"),
                         n_trials_per_cond = 2,
                         n_item_class = 2, n_item_per_class = 5,
                         item_class_weights, sd_bA = 0.2,
                         b_stick = b_stick, sd_b_stick = 1,
                         rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                         rho_psi = rho_psi, sd_rho_psi = 0.25,
                         abs_dir_tuning = abs_dir_tuning,
                         inital_sel_params = inital_sel_params,
                         rel_proximity = FALSE,
                         filename = "test_anna") 

m <- fit_model(d, fomo_ver = "1.0", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.1", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.2", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.3", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.4", mode = "traintest", iter = 500) 
m <- fit_model(d, fomo_ver = "1.5", mode = "traintest", iter = 500) 
