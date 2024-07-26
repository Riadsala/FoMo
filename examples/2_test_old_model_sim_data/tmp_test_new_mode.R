
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

item_class_weights = list(c(0.7, 0.3, 0, 0))

b_stick = 1

rho_delta = 15
sd_rho_delta = 5

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

d2 <- sim_foraging_people(n_people = 5,
                          n_conditions = 1,
                          cond_lab = c("simple test"),
                          n_trials_per_cond = 4,
                          n_item_class = 2, n_item_per_class = 20,
                          item_class_weights, sd_bA = 0.2,
                          b_stick = b_stick, sd_b_stick = 1,
                          rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                          rho_psi = rho_psi, sd_rho_psi = 0.5,
                          abs_dir_tuning = abs_dir_tuning,
                          inital_sel_params = inital_sel_params) 

d2$found <- fix_person_and_trial(d2$found)
d2$stim <- fix_person_and_trial(d2$stim)


d2_list <- prep_data_for_stan(d2$found, d2$stim, c("spatial", "item_class"))
d2_list_p <- add_priors_to_d_list(d2_list, modelver = "1.0")
d2_list_p$n_trials_to_sim <- 10


iter = 500
mod <- cmdstan_model("../../models/multi_level/FoMo1_0_2.stan", 
                     cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)

fit <- mod$sample(data = d2_list_p, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 200, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)
