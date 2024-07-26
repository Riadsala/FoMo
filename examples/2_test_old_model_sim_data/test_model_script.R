library(tidyverse)
library(cmdstanr)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")

################################################
# First test single level model

n_trials_per_cond <- 20

n_item_class <- 2
n_item_per_class <- 20
item_class_weights = c(0.5, 0.5, 0, 0)
b_stick = 2
b_memory = 0

abs_dir_tuning = list(kappa = rep(20, 4), theta = c(2, 0.5, 1, 0.5))
rho_delta = 15
rho_psi = 5

d <- sim_foraging_multiple_trials(person = 1, 
                                   condition = "test",
                                   n_trials_per_cond = n_trials_per_cond,
                                   n_item_class =  n_item_class, n_item_per_class = n_item_per_class,
                                   item_class_weights = item_class_weights, item_labels = item_labels,
                                   b_stick = b_stick, 
                                   rho_delta = rho_delta, 
                                   rho_psi = rho_psi, 
                                   abs_dir_tuning = abs_dir_tuning,
                                   b_memory = b_memory,
                                   inital_sel_params = inital_sel_params,
                                   init_sel_lambda = init_sel_lambda)

d_list <- prep_data_for_stan(d$found, d$stim, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.0")
d_list$n_trials_to_sim <- 10

iter = 500
mod <- cmdstan_model("../../models/simple/FoMo1_0.stan", 
                     cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 200, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

post <- extract_post(fit, d, multi_level = FALSE)

plot_model_fixed(post,  gt = list(b_a = qlogis(item_class_weights[1]),
                                  b_stick = b_stick,
                                  rho_delta = rho_delta,
                                  rho_psi = rho_psi))


pred <- summarise_postpred(fit, d, multi_level = FALSE, get_sim = FALSE)

pred$acc %>% group_by(found, .draw) %>%
  summarise(model_correct = mean(model_correct)) %>%
  summarise(model_correct = mean(model_correct)) %>%
  ggplot(aes(found, model_correct)) + 
  geom_path()
  


################################################
# Now test multi-level model

d <- sim_foraging_people(n_people = 5,
                          n_conditions = 2,
                          cond_lab = c("anna", "banana"),
                          n_trials_per_cond = 10,
                          n_item_class = 2, n_item_per_class = 10,
                          list(item_class_weights, item_class_weights), sd_bA = 0.1,
                          b_stick = c(0, 1), sd_b_stick = 0.1,
                          rho_delta = c(15, 15), sd_rho_delta = 1,
                          rho_psi = c(0, 0), sd_rho_psi = 0.1,
                          abs_dir_tuning = abs_dir_tuning,
                          inital_sel_params = inital_sel_params) 

d$found <- fix_person_and_trial(d$found)
d$stim <- fix_person_and_trial(d$stim)

d_list <- prep_data_for_stan(d$found, d$stim, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.0")
d_list$n_trials_to_sim <- 10

iter = 500
mod <- cmdstan_model("../../models/multi_level/FoMo1_0.stan", 
                     cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 200, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

post <- extract_post(fit, d, multi_level = TRUE)
plot_model_fixed(post,   gt = list(b_a = qlogis(item_class_weights[1]),
                                   b_stick = b_stick,
                                   rho_delta = rho_delta,
                                   rho_psi = rho_psi))


pred <- summarise_postpred(fit, d, multi_level = TRUE, get_sim = TRUE)

pred$acc %>% group_by(found, .draw) %>%
  summarise(model_correct = mean(model_correct)) %>%
  summarise(model_correct = mean(model_correct)) %>%
  ggplot(aes(found, model_correct)) + 
  geom_path()
