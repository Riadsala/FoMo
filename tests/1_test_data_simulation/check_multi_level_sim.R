library(tidyverse)

source("../../functions/sim_foraging_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/prep_data.R")

######################################################################
# lets simulate some data 
######################################################################

experiment_params <- list(n_people = 16, 
                          n_conditions = 2,
                          condition_labels = c("A", "B"),
                          n_trials_per_cond = 2)

stimuli_params <- list(n_item_class = 4,
                       n_item_per_class = c(10, 10, 10, 10), 
                       item_labels = c("a", "b", "d1", "d2"))

foraging_params <- list(b_a = c(-1, 1), 
                        b_s = c(1, -1), 
                        rho_delta = c(0.5, 1.5),
                        rho_psi = c(-0.5, 0.5))

variance_params <- list(b_a = c(0.01, 0.01), 
                        b_s = c(0.01, 0.01), 
                        rho_delta = c(0.01, 0.01),
                        rho_psi = c(0.01, 0.01))

absdir_params <- "off"
initsel_params <- "off"

params <- list(e = experiment_params,
               s = stimuli_params,
               f = foraging_params,
               v = variance_params,
               a = absdir_params,
               i = initsel_params)

d <- sim_foraging_people(params) 

######################################################################
# now check things worked as expected
######################################################################

# first, plot the individual level parameters
foraging_params %>% as_tibble() %>%
  mutate(condition = experiment_params$condition_labels) %>%
  pivot_longer(-condition, names_to = "param") -> gt

d$dp %>% pivot_longer(-c(person, condition), names_to = "param") %>%
  ggplot(aes(value, fill = condition)) + 
  geom_histogram(alpha = 0.5) +
  geom_vline(data = gt, aes(xintercept = value, colour = condition)) + 
  facet_wrap(~param, scales = "free")

# now check run statistics
rl <- get_run_info_over_trials(d$found) 

rl %>% select(person, condition, max_run_length, n_runs) %>%
  pivot_longer(-c(person, condition), names_to = "stat") %>%
  ggplot(aes(value, fill = condition)) + 
  geom_histogram(alpha = 0.5, binwidth = 1, position = position_identity()) + 
  facet_wrap(~ stat, scales = "free")

# check sim param correlates with run statistics
rl %>% group_by(person, condition) %>%
  summarise(max_run_length = median(max_run_length),
            n_runs = median(n_runs)) %>%
  full_join(d$dp) %>%
  pivot_longer(c("b_a", "b_s", "rho_delta", "rho_psi"), names_to = "param", values_to = "x") %>%
  pivot_longer(c("max_run_length", "n_runs"), names_to = "stat", values_to = "y")  %>%
  ggplot(aes(x, y, colour = condition)) +
  geom_point() + 
  facet_grid(param~stat)

######################################################################
# what about prep_data() ? 
######################################################################

dl <- prep_data_for_stan(d)

head(t(rbind(dl$Z, dl$X)), 20)
dl$trial
