library(tidyverse)
library(cmdstanr)
library(patchwork)

source("../functions/import_data.R")
source("../functions/prep_data.R")
source("../functions/compute_summary_stats.R")
source("../functions/plot_model.R")
source("../functions/plot_data.R")
source("../functions/post_functions.R")
source("../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

model_ver <- "1_5"
dataset <- "hughes2024rsos"

# read in data
d <- import_data(dataset)

folder <- paste0("1_fit_models/scratch/post/", dataset, "/")

plot_models_accuracy(dataset)

#############################################################################
# plot accuracy
#############################################################################
acc <- read_csv(paste0(folder, "acc_train1_0.csv"))
acc_plt <- plot_model_accuracy(acc)
rm(acc)

#############################################################################
# plot posterior densities
#############################################################################
m <- readRDS(paste0("1_fit_models/scratch/models/", dataset, "/train", model_ver, ".model"))
post <- extract_post(m, d)
post_plt <- plot_model_fixed(post)

#############################################################################
# create plot
#############################################################################
acc_plt / post_plt

#############################################################################
# compute & compare run statistics
#############################################################################
pred <- readRDS(paste0(folder, "pred_train1_0.rds"))

rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs       = mean(n_runs),
            .groups = "drop")

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$trialwise %>% filter(.draw == 1)) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs       = mean(n_runs))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) %>%
  pivot_longer(c(max_run_length, num_runs), names_to = "statistic") %>%
  pivot_wider(names_from = "x") -> rl

# plot
rl %>% ggplot(aes(predicted, observed, 
                  colour = condition, shape = condition)) + 
  geom_point() + 
  geom_abline(linetype = 2) + 
  facet_wrap(~statistic, scales = "free")

#############################################################################
# compute & compare iisv statistics
#############################################################################

rle <- get_iisv_over_trials(d$found) 

# compute simulated run statistics
rlp <- get_iisv_over_trials(pred$trialwise %>% filter(.draw == 1))

  # bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) -> rl

# plot
rl %>% ggplot(aes(predicted, observed, 
                  colour = condition, shape = condition)) + 
  geom_point() + 
  geom_abline(linetype = 2) + 
  facet_wrap(~statistic, scales = "free")
