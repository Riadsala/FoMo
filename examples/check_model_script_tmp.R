library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

d <- import_data("kristjansson2014plos")
m <- readRDS("scratch/kristjansson2014plos_train_1_0.model")
post <- extract_post(m, d)
plot_model_fixed(post)


# get simulation data for model
pred <- summarise_postpred(list(training = m, testing = m), d, 
                           get_sim = TRUE, draw_sample_frac = 0.001)

pred$acc %>% group_by(.draw, condition, found) %>%
  summarise(correct = mean(model_correct)) %>%
  ggplot(aes(found, correct, fill = condition)) + 
  stat_lineribbon(alpha = 0.3)

# compute empirical run statistics
rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$sim) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) %>%
  pivot_wider(names_from = "x", values_from = "max_run_length") -> rl

ggplot(rl, aes(observed, predicted, colour = condition)) + geom_point() +
  geom_abline()


##########################################

# compute empirical run statistics
rle <- get_iisv_over_trials(d$found) %>%
  filter(is.finite(d2)) %>%
  group_by(condition, found) %>%
  summarise(d2 = mean(d2))

# compute simulated run statistics
rlp <- get_iisv_over_trials(pred$sim) %>%
  filter(is.finite(d2)) %>%
  group_by(condition, found) %>%
 summarise(d2 = mean(d2))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) -> iisv

ggplot(iisv, aes(found, d2, colour = x)) + geom_point() +
  geom_abline()  +
  facet_wrap(~condition) + 
  scale_y_log10()



