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

model_ver <- "1_0"
# dataset   <- "kristjansson2014plos"
dataset <- "clarke2022qjep"

# read in data
d <- import_data(dataset)

# read in model and predictions for test data
m <- read_rds(paste0("scratch/", dataset, "_train_", model_ver, ".model"))
t <- read_rds(paste0("scratch/", dataset, "_test_", model_ver, ".model"))

# get simulation data for model
pred <- summarise_postpred(list(training = m, testing = t), d, 
                           get_sim = TRUE, draw_sample_frac = 0.01)

# sanity check we have correct number of items per trial etc

pred$sim %>% group_by(person, condition, .draw, trial) %>%
  summarise(n = n())  -> tmp


pred$acc %>% group_by(.draw) %>%
  summarise(n = n()) %>%
  summarise(n = unique(n))


 tmp %>% filter(n > 40) -> tmp2
 
 summary(tmp2)
 
 filter(pred$acc, .draw == 1397) %>% group_by(person, trial, condition) %>%
   summarise(n=n()) -> tmp3

rm(m, t)



#############################################################################

# compute empirical run statistics
rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs       = mean(n_runs))

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$sim %>% filter(.draw == 1101)) %>%
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

