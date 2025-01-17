library(tidyverse)
library(cmdstanr)
library(loo)

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

item_class_weights = list(c(0.5, 0.5, 0, 0),
                          c(0.7, 0.3, 0, 0))

b_stick = c(0, 2)

rho_delta = c(20, 15)
sd_rho_delta = 5

rho_psi = c(-1, -1)

abs_dir_tuning = list(kappa = rep(10, 4), theta = rep(1, 4))

d <- sim_foraging_people(n_people = 8,
                          n_conditions = 2,
                          cond_lab = c("A", "B"),
                          n_trials_per_cond = 4,
                          n_item_class = 2, n_item_per_class = 9,
                          item_class_weights, sd_bA = 0.2,
                          b_stick = b_stick, sd_b_stick = 1,
                          rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                          rho_psi = rho_psi, sd_rho_psi = 0.5,
                          abs_dir_tuning = abs_dir_tuning,
                          inital_sel_params = inital_sel_params)

d$found <- fix_person_and_trial(d$found)
d$stim <- fix_person_and_trial(d$stim)

d_list <- prep_data_for_stan(d, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.0")

iter = 100
mod <- cmdstan_model("../../models/multi_level/FoMo1_0.stan")

fit <- mod$sample(data = d_list, 
                  chains = 1, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

# compute run statistics for later
runs_emp <- get_run_info_over_trials(d$found) 
iisv_emp <- get_iisv_over_trials(d$found) 

post <- extract_post(fit, d, multi_level = TRUE)
plot_model_fixed(post)

pred <- summarise_postpred(fit, d, draw_sample_frac = 0.1)

pred$sim %>% group_by(.draw, person, trial) %>%
  summarise(n = n()) %>%
  summary

###############################################

simple_run_stat_comparison <- function(runs_emp, pred) {
  
  runs_sim <- get_run_info_over_trials(pred$sim) 
  
  bind_rows(runs_emp %>% mutate(class = "empirical"), 
            runs_sim %>% mutate(class = "simulated") %>% select(-.draw)) %>%
    pivot_longer(c("max_run_length", "n_runs")) %>%
    ggplot(aes(class, value, fill = condition)) + 
    geom_boxplot() + 
    facet_wrap(~name, scales = "free")
}

simple_run_stat_comparison(runs_emp, pred)
