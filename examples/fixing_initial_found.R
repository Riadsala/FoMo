library(tidyverse)
library(cmdstanr)

# does fixining the first item to be found improve path metrics?

source("../functions/import_data.R")
source("../functions/compute_summary_stats.R")
source("../functions/post_functions.R")

dataset <- "clarke2022qjep"

d <- import_data(dataset)
folder <- paste0("1_fit_models/scratch/models/", dataset, "/sim/")

########################################################
# compute postpred stats for fif version of the genquant
########################################################

# get fif model predictions

modelver <- "1_3fif"
predfif <- extract_pred(dataset, modelver, folder)

# get run info for fif data
rl_fif <- get_run_info_over_trials(predfif$trialwise %>% filter(.draw == 1)) %>%
  group_by(.draw, person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs = mean(n_runs),
            mean_bestr = mean(best_r),
            mean_pao = mean(pao),
            .groups = "drop") %>% 
  pivot_longer(-c(.draw, person, condition), names_to = "statistic", values_to = "fif")

# now join with standard model
rl <- read_csv("1_fit_models/scratch/post/tagu2022cog/run_statistics.csv")

rl %>% select(-v1_0, -v1_2) %>%
  full_join(rl_fif) %>%
  pivot_longer(c(v1_3, fif), names_to = "init_sel") -> rl

rl %>% ggplot(aes(observed, value, colour = condition)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_abline() + 
  ggh4x::facet_grid2(init_sel ~ statistic, scale = "free", independent = "all")

rl %>% mutate(abs_err = abs(value - observed)) %>%
  group_by(condition, statistic, init_sel) %>%
  summarise(mean_abserr = mean(abs_err)) %>%
  pivot_wider(names_from = init_sel, values_from = "mean_abserr") %>%
  mutate(improvement = fif-v1_3,
         prop = (v1_3-fif)/v1_3)
