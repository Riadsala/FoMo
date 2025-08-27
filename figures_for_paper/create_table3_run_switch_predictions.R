library(tidyverse)
library(tidybayes)
source("../functions/import_data.R")
# script to compute switch acc

dataset <- "kristjansson2014plos"
model_ver <- "1_0"

p <- readRDS(paste0("../examples/1_fit_models/scratch/post/", dataset, "/pred_", model_ver, ".rds"))

p <- p$itemwise %>%
  select(-model_correct, -x, -y, -trial)

# add in p_item_class
d <- import_data(dataset)
d$stim %>%
  select(-trial_p) %>%
  select(person, condition, trial_p = "trial", P = "id", pred_class = "item_class") %>%
  right_join(p) %>%
  select(-P, -id) %>% 
  arrange(.draw, person, condition, trial_p, found) %>%
  select(.draw, person, condition, trial_p, found, item_class, pred_class) %>%
  mutate(
    prev_class = lag(item_class),
    e_switch = item_class != prev_class,
    p_switch = pred_class != prev_class,
    pred_correct = e_switch == p_switch) %>%
  filter(found > 1) -> p


p %>% group_by(condition, e_switch,.draw, person, trial_p) %>%
  summarise(mean_switch_acc = mean(pred_correct)) %>%
  summarise(mean_switch_acc = mean(mean_switch_acc)) %>%
  summarise(mean_switch_acc = mean(mean_switch_acc)) %>%
  median_hdci(mean_switch_acc)
