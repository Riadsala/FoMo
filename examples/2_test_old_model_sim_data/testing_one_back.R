library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")

options(mc.cores = 1)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())


n_trials_per_cond <- 25

n_item_class <- 2
n_item_per_class <- 20
item_class_weights = c(0.7, 0.3, 0, 0)
b_stick = 2
b_memory = 0

abs_dir_tuning = list(kappa = rep(20, 4), theta = c(0, 00, 0, 0))
rho_delta = 10

dbm <- tibble()

for (rho_psi in seq(-3, 2, 0.5)) {

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
  
  
  iter = 500
  mod <- cmdstan_model("../../models/simple/FoMo1_x.stan", 
                       cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)
  
  # run model
  m <- mod$sample(data = d1_list, 
                             chains = 4, parallel_chains = 4, threads = 4,
                             refresh = 0, 
                             iter_warmup = iter, iter_sampling = iter,
                             sig_figs = 3)
  
  post <- extract_post(m ,d1, multi_level = FALSE)
  
  post$fixed %>% select(.draw, b_m) %>%
    mutate(rho_psi = rho_psi) -> post
  
  dbm %>% bind_rows(post) -> dbm

}


dbm %>% group_by(rho_psi) %>%
  median_hdci(exp(b_m), .width = c(0.53, 0.97)) %>%
  ggplot(aes(x = rho_psi, ymin = .lower, ymax = .upper, group = .width)) +
  geom_ribbon(alpha = 0.5) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) 
  


############ Look at some people from Tagu?


d <- import_data('tagu2022cog')

mod11 <- cmdstan_model("../../models/simple/FoMo1_1.stan", 
                       cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)


iter = 500


pp = 6


  d_one_person <- filter_one_person(d, pp) 
  
  d_list <- prep_data_for_stan(d_one_person$found, d_one_person$stim, c("spatial", "item_class"))
  # add priors to list
  d_list <- add_priors_to_d_list(d_list)
  d_list$n_trials_to_sim <- 1
  
  filename <- paste0("scratch/person", pp)
  
  # model 1.1
  fit <- mod11$sample(data = d_list, 
                      chains = 4, parallel_chains = 4, threads = 4,
                      refresh = 100, 
                      init = 1,
                      iter_warmup = iter, iter_sampling = iter,
                      sig_figs = 3)
  
 

 
  # model 1.1
  fitX <- mod$sample(data = d_list, 
                      chains = 4, parallel_chains = 4, threads = 4,
                      refresh = 100, 
                      init = 1,
                      iter_warmup = iter, iter_sampling = iter,
                      sig_figs = 3)
  
  
  loo1 <- fit$loo()
  looX <- fitX$loo()  

  
  loo::loo_model_weights(list("FoMo 1.1" = loo1, "FoMo 1.X" = looX))  
  