library(tidyverse)
library(cmdstanr)
library(loo)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")

# create scratch folder if it does not yet exist
if (!file.exists("scratch"))  dir.create("scratch")

## model 1.0

item_class_weights = list(c(0.7, 0.3, 0, 0))

b_stick = 1

rho_delta = 20
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

d <- sim_foraging_people(n_people = 16,
                         n_conditions = 1,
                         cond_lab = c("simple test"),
                         n_trials_per_cond = 8,
                         n_item_class = 2, n_item_per_class = 10,
                         item_class_weights, sd_bA = 0.2,
                         b_stick = b_stick, sd_b_stick = 1,
                         rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                         rho_psi = rho_psi, sd_rho_psi = 0.25,
                         abs_dir_tuning = abs_dir_tuning,
                         inital_sel_params = inital_sel_params) 

d$found <- fix_person_and_trial(d$found)
d$stim <- fix_person_and_trial(d$stim)

saveRDS(d, "scratch/d_1_0.rds")


m <- fit_model(d, fomo_ver = "1.0", mode = "all",  iter = 500) 

  ## model 1.1

item_class_weights = list(c(0.7, 0.3, 0, 0))

b_stick = 1

rho_delta = 1
sd_rho_delta = 0.1

rho_psi = -1

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

d2 <- sim_foraging_people(n_people = 10,
                          n_conditions = 1,
                          cond_lab = c("simple test"),
                          n_trials_per_cond = 4,
                          n_item_class = 2, n_item_per_class = 10,
                          item_class_weights, sd_bA = 0.2,
                          b_stick = b_stick, sd_b_stick = 1,
                          rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                          rho_psi = rho_psi, sd_rho_psi = 0.5,
                          abs_dir_tuning = abs_dir_tuning,
                          inital_sel_params = inital_sel_params,
                          rel_proximity = TRUE) 

d2$found <- fix_person_and_trial(d2$found)
d2$stim <- fix_person_and_trial(d2$stim)

saveRDS(d2, "scratch/d_1_1.rds")

iter = 200
mod <- cmdstan_model("../../models/multi_level/FoMo1_1.stan")

d2_list <- prep_data_for_stan(d2$found, d2$stim, c("spatial", "item_class"))
d2_list <- add_priors_to_d_list(d2_list, modelver = "1.1")
d2_list$n_trials_to_sim <- 1

fit <- mod$sample(data = d2_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object("scratch/multi_level_1_1_tmp.rds")

# model 1.2

mod <- cmdstan_model("../../models/multi_level/FoMo1_2.stan")

fit <- mod$sample(data = d2_list, 
                           chains = 4, parallel_chains = 4, threads = 4,
                           refresh = 10, 
                           iter_warmup = iter, iter_sampling = iter,
                           sig_figs = 3)

fit$save_object("scratch/multi_level_1_2_tmp.rds")


### model testing using train-test split

# going to use d2

# model 1.0

df_train <- d2$found %>%
  filter((trial %% 2) == 0) %>%
  group_by(trial) %>%
  mutate(trial = cur_group_id()) %>%
  ungroup()

ds_train <- d2$stim %>%
  filter((trial %% 2) == 0) %>%
  group_by(trial) %>%
  mutate(trial = cur_group_id()) %>%
  ungroup()

d_train_list <- prep_data_for_stan(df_train, ds_train, c("spatial", "item_class"))
d_train_list <- add_priors_to_d_list(d_train_list, modelver = "1.0")
d_train_list$n_trials_to_sim <- 1

iter = 100
mod <- cmdstan_model("../../models/multi_level/FoMo1_0.stan")

# run model
m_train_1_0 <- mod$sample(data = d_train_list, 
                          chains = 4, parallel_chains = 4, threads = 4,
                          refresh = 0, 
                          iter_warmup = iter, iter_sampling = iter,
                          sig_figs = 3)

df_test <- d2$found %>%
  filter((trial %% 2) != 0) %>%
  group_by(trial) %>%
  mutate(trial = cur_group_id()) %>%
  ungroup()

ds_test <- d2$stim %>%
  filter((trial %% 2) != 0) %>%
  group_by(trial) %>%
  mutate(trial = cur_group_id()) %>%
  ungroup()

d_test_list <- prep_data_for_stan(df_test, ds_test, model_components = c("spatial", "item_class"))
d_test_list <- add_priors_to_d_list(d_test_list, modelver = "1.0")
d_test_list$n_trials_to_sim <- 1


gen_test <- mod$generate_quantities(m_train_1_0, data = d_test_list, seed = 123)
log_pd_kfold <- gen_test$draws("log_lik", format = "matrix")

elpd_kfold_1_0 <- elpd(log_pd_kfold)
saveRDS(elpd_kfold_1_0, "scratch/elpd_1_0.rds")

## model 1.1

d_train_list <- add_priors_to_d_list(d_train_list, modelver = "1.1")

mod <- cmdstan_model("../../models/multi_level/FoMo1_1.stan")

# run model
m_train_1_1 <- mod$sample(data = d_train_list, 
                          chains = 4, parallel_chains = 4, threads = 4,
                          refresh = 10, 
                          iter_warmup = iter, iter_sampling = iter,
                          sig_figs = 3)

d_test_list <- add_priors_to_d_list(d_test_list, modelver = "1.1")

gen_test <- mod$generate_quantities(m_train_1_1, data = d_test_list, seed = 123)
log_pd_kfold <- gen_test$draws("log_lik", format = "matrix")

elpd_kfold_1_1 <- elpd(log_pd_kfold)
saveRDS(elpd_kfold_1_1, "scratch/elpd_1_1.rds")

## model 1.2

mod <- cmdstan_model("../../models/multi_level/FoMo1_2.stan")

# run model
m_train_1_2 <- mod$sample(data = d_train_list, 
                          chains = 4, parallel_chains = 4, threads = 4,
                          refresh = 10, 
                          iter_warmup = iter, iter_sampling = iter,
                          sig_figs = 3)

gen_test <- mod$generate_quantities(m_train_1_2, data = d_test_list, seed = 123)
log_pd_kfold <- gen_test$draws("log_lik", format = "matrix")

elpd_kfold_1_2 <- elpd(log_pd_kfold)
saveRDS(elpd_kfold_1_2, "scratch/elpd_1_2.rds")


#### Simulating multiple conditions

item_class_weights = list(c(0.5, 0.5, 0, 0), 
                          c(0.7, 0.3, 0, 0))

b_stick = c(0, 2)

rho_delta = c(20, 15)
sd_rho_delta = 5

rho_psi = c(-1, -1)

abs_dir_tuning = list(kappa = rep(10, 4), theta = rep(1, 4))

d3 <- sim_foraging_people(n_people = 10,
                         n_conditions = 2,
                         cond_lab = c("A", "B"),
                         n_trials_per_cond = 4,
                         n_item_class = 2, n_item_per_class = 10,
                         item_class_weights, sd_bA = 0.2,
                         b_stick = b_stick, sd_b_stick = 1,
                         rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                         rho_psi = rho_psi, sd_rho_psi = 0.5,
                         abs_dir_tuning = abs_dir_tuning,
                         inital_sel_params = inital_sel_params) 

d3$found <- fix_person_and_trial(d3$found)
d3$stim <- fix_person_and_trial(d3$stim)

saveRDS(d3, "scratch/d_2cond.rds")

d_list <- prep_data_for_stan(d3$found, d3$stim, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.0")
d_list$n_trials_to_sim <- 1

iter = 200
mod <- cmdstan_model("../../models/multi_level/FoMo1_0.stan")

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object("scratch/multi_level_2cond_1_0_tmp.rds")

## model 1.1

rho_delta = c(1,2)
sd_rho_delta = 0.1

d4 <- sim_foraging_people(n_people = 10,
                          n_conditions = 2,
                          cond_lab = c("A", "B"),
                          n_trials_per_cond = 4,
                          n_item_class = 2, n_item_per_class = 10,
                          item_class_weights, sd_bA = 0.2,
                          b_stick = b_stick, sd_b_stick = 1,
                          rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                          rho_psi = rho_psi, sd_rho_psi = 0.5,
                          abs_dir_tuning = abs_dir_tuning,
                          inital_sel_params = inital_sel_params,
                          rel_proximity = TRUE) 

d4$found <- fix_person_and_trial(d4$found)
d4$stim <- fix_person_and_trial(d4$stim)

saveRDS(d4, "scratch/d_2cond_relprox.rds")

d_list <- prep_data_for_stan(d4$found, d4$stim, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.1")
d_list$n_trials_to_sim <- 1

iter = 200
mod <- cmdstan_model("../../models/multi_level/FoMo1_1.stan")

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object("scratch/multi_level_2cond_1_1_tmp.rds")

# model 1.2

mod <- cmdstan_model("../../models/multi_level/FoMo1_2.stan")

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object("scratch/multi_level_2cond_1_2_tmp.rds")