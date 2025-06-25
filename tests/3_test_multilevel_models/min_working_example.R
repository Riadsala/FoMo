library(tidyverse)
library(cmdstanr)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")

options(mc.cores = 8)

######################################################################
# lets simulate some data - this is for the TEST MULTICOND dataset
######################################################################

experiment_params <- list(n_people = 10, 
                          n_conditions = 2,
                          condition_labels = c("A", "B"),
                          n_trials_per_cond = 6)

stimuli_params <- list(n_item_class = 4,
                       n_item_per_class = c(20, 20, 10, 10), 
                       item_labels = c("a", "b", "d1", "d2"))

foraging_params <- list(b_a = c(0, 0.25), 
                        b_s = c(0.5, 2.5), 
                        rho_delta = c(1, 0.8),
                        rho_psi = c(-1, 1))

variance_params <- list(b_a = c(0.1, 0.1), 
                        b_s = c(0.1, 0.1), 
                        rho_delta = c(0.1, 0.1),
                        rho_psi = c(0.1, 0.1))

absdir_params <- list(
  kappa = rep(20, 4), theta = c(5, 1, 5, 1))

# absdir_params <- "off"

initsel_params <- "off"

params <- list(e = experiment_params,
               s = stimuli_params,
               f = foraging_params,
               v = variance_params,
               a = absdir_params,
               i = initsel_params)

d <- sim_foraging_people(params) 

plot_a_trial(d$stim, d$found, 1)

######################################################################
# fit model
######################################################################

dl <- prep_data_for_stan(d)
dl <- add_priors_to_d_list(dl, modelver = "1.3")


iter = 500
mod <- cmdstan_model("../../models/multi_level/FoMo1_3.stan",
                     cpp_options = list(stan_threads = TRUE))

m <- mod$sample(data = dl, 
                iter_warmup = iter, iter_sampling = iter,
                chains = 4, 
                parallel_chains = 8,
                threads_per_chain = 2)

######################################################################
# check model posterior
######################################################################

m$summary()

post <- extract_post(m, d)
plot_model_fixed(post, gt = params)

plot_model_random(post)
# 
# plot_model_theta(post)

bayesplot::mcmc_trace(m$draws(), par = "lp__")
