library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")

options(mc.cores = 4)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())


mod12 <- cmdstan_model("../../models/simple/FoMo1_2.stan", 
                       cpp_options = list(stan_threads = TRUE))

iter = 500

n_trials_per_cond <- 25

n_item_class <- 2
n_item_per_class <- 20
item_class_weights = c(0.7, 0.3, 0, 0)
b_stick = 1
b_memory = 0

abs_dir_tuning = list(kappa = rep(20, 4), theta = c(0, 0, 0, 0))
rho_delta = 5

bms <- tibble()


for (rho_psi in seq(-3, 2, 2)) {
  
  
  d <- sim_foraging_multiple_trials(person = 1, 
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
  
  d_list <- prep_data_for_stan(d$found, d$stim, c("spatial", "item_class"))
  d_list <- add_priors_to_d_list(d_list, modelver = "1.1")
  d_list$n_trials_to_sim <- 1
  
  # model 1.1
  fit <- mod12$sample(data = d_list, 
                      chains = 4, parallel_chains = 4, threads = 4,
                      refresh = 100, 
                      init = 1,
                      iter_warmup = iter, iter_sampling = iter,
                      sig_figs = 3)
  
  p <- extract_post(fit, d, multi_level = FALSE)
  
  p$fixed$b_m %>% 
    median_hdci(b_m) %>% 
    mutate(rho_psi = rho_psi) %>%
    as_tibble() %>%
    add_row(bms) -> bms
}


bms %>% ggplot(aes(rho_psi, exp(y))) + 
  geom_path() +
  geom_pointinterval(aes(ymin = exp(ymin), ymax = exp(ymax)))
