library(tidyverse)
library(cmdstanr)
library(patchwork)

source("../../functions/import_data.R")
source("../../functions/post_functions.R")
source("../../functions/plot_model.R")

d <- import_data("tagu2025")

m <- readRDS("../1_fit_models/scratch/models/tagu2025all1_0.model")

# First of all, plot accuracy
pred <- extract_pred(m, d)
acc <- summarise_acc(pred)
plot_model_accuracy(acc)


# Posterior Density Distributions
post <- extract_post(m, d)
plot_model_fixed(post)


# are there differences between conditions?
post$fixed %>% 
  pivot_longer(-c(.draw, condition), names_to = "param") %>%
  pivot_wider(names_from = "condition") %>%
  group_by(.draw, param) %>%
  mutate(negative = negative - neutral,
         positive = positive - neutral,
         .keep = "none") %>%
  pivot_longer(c(negative, positive), names_to = "comparison") %>%
  ggplot(aes(value, fill = comparison)) + 
  geom_density(alpha = 0.6) + 
  geom_vline(xintercept = 0, linetype = 2) + 
  facet_wrap(~param, scales = "free")
  
# Look at individuals 
plot_model_random(post)
