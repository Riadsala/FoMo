#!/usr/bin/Rscript

library(cmdstanr)
library(tidyverse)

options(mc.cores = 1)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")

#########################
# now fit models
########################

d <- import_data('tagu2022cog')

d_list <- prep_data_for_stan(d$found, d$stim, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.0")

iter = 1000
d_list$n_trials_to_sim <- 1

mod <- cmdstan_model("../../models/multi_level/FoMo1_0.stan", 
                     cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  init = 1,
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object("scratch/tagu_1_0_tmp.rds")

# pre save moment match loo
fit1_0_loo <- fit$loo(moment_match = TRUE)
saveRDS(fit1_0_loo, "scratch/tagu_1_0_loo.rds")