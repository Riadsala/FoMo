library(tidyverse)
library(cmdstanr)

options(mc.cores = 8, 
        digits = 2)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())


source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")


mod_old <- cmdstan_model("../../models/simple/FoMo1_0.stan", 
                     cpp_options = list(stan_threads = TRUE))

mod_new <- cmdstan_model("../../models/simple/FoMo1_0_2.stan", 
                     cpp_options = list(stan_threads = TRUE))





n_trials_per_cond <- 10

n_item_class <- 2
n_item_per_class <- 20
item_class_weights = c(0.7, 0.3, 0, 0)
b_stick = 2
b_memory = 0

abs_dir_tuning = list(kappa = rep(20, 4), theta = c(2, 0.5, 1, 0.5))
rho_delta = 10
rho_psi = 5

d1 <- sim_foraging_multiple_trials(person = 1, 
                                   condition = "test",
                                   n_item_class =  n_item_class, n_item_per_class = n_item_per_class,
                                   item_class_weights = item_class_weights, item_labels = item_labels,
                                   b_stick = b_stick, 
                                   rho_delta = rho_delta, 
                                   rho_psi = rho_psi, 
                                   abs_dir_tuning = abs_dir_tuning,
                                   b_memory = b_memory,
                                   inital_sel_params = inital_sel_params,
                                   init_sel_lambda = init_sel_lambda)



d1_list <- prep_data_for_stan(d1$found, d1$stim, c("spatial", "item_class"))

# add priors to list
d1_list$prior_mu_b_a <- 0
d1_list$prior_sd_b_a <- 0.5
d1_list$prior_mu_b_stick <- 0
d1_list$prior_sd_b_stick <- 1
d1_list$prior_mu_rho_delta <- 15
d1_list$prior_sd_rho_delta <- 5
d1_list$prior_mu_rho_psi <- 0
d1_list$prior_sd_rho_psi <- 1
d1_list$n_trials_to_sim <- 10



iter = 1000
# run model
m_simple_1_old <- mod_old$sample(data = d1_list, 
                           chains = 4, parallel_chains = 4, threads = 4,
                           refresh = 500, 
                           iter_warmup = iter, iter_sampling = iter,
                           sig_figs = 3)

m_simple_1_new <- mod_new$sample(data = d1_list, 
                                 chains = 4, parallel_chains = 4, threads = 4,
                                 refresh = 500, 
                                 iter_warmup = iter, iter_sampling = iter,
                                 sig_figs = 3)


post_old <- extract_post(m_simple_1_old, d1, multi_level = FALSE)
post_new <- extract_post(m_simple_1_new, d1, multi_level = FALSE)



# plot model
plot_model_fixed(post_old, gt = list(b_a = plogis(item_class_weights[1]),
                                 b_stick = b_stick,
                                 rho_delta = rho_delta,
                                 rho_psi = rho_psi))


# plot model
plot_model_fixed(post_new, gt = list(b_a = plogis(item_class_weights[1]),
                                     b_stick = b_stick,
                                     rho_delta = rho_delta,
                                     rho_psi = rho_psi))
