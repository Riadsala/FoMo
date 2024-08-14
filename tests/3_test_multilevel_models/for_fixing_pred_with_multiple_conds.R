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

# this is how the data were generated

# item_class_weights = list(c(0.5, 0.5, 0, 0), 
#                           c(0.7, 0.3, 0, 0))
# 
# b_stick = c(0, 2)
# 
# rho_delta = c(20, 15)
# sd_rho_delta = 5
# 
# rho_psi = c(-1, -1)
# 
# abs_dir_tuning = list(kappa = rep(10, 4), theta = rep(1, 4))
# 
# d3 <- sim_foraging_people(n_people = 10,
#                           n_conditions = 2,
#                           cond_lab = c("A", "B"),
#                           n_trials_per_cond = 4,
#                           n_item_class = 2, n_item_per_class = 10,
#                           item_class_weights, sd_bA = 0.2,
#                           b_stick = b_stick, sd_b_stick = 1,
#                           rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
#                           rho_psi = rho_psi, sd_rho_psi = 0.5,
#                           abs_dir_tuning = abs_dir_tuning,
#                           inital_sel_params = inital_sel_params) 
# 
# d3$found <- fix_person_and_trial(d3$found)
# d3$stim <- fix_person_and_trial(d3$stim)

d <- readRDS("scratch/d_2cond.rds")
fit <- readRDS("scratch/multi_level_2cond_1_0_tmp.rds")

item_class_weights = list(c(0.5, 0.5, 0, 0), 
                          c(0.7, 0.3, 0, 0))

b_stick = c(0, 2)
rho_delta = c(20, 15)
rho_psi = c(-1, -1)

# compute run statistics for later
runs_emp <- get_run_info_over_trials(d$found) 
iisv_emp <- get_iisv_over_trials(d$found) 

post <- extract_post(fit, d, multi_level = TRUE)

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