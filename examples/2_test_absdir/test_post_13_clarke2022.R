library(tidyverse)
library(cmdstanr)
library(tidybayes)

source("../../functions/post_functions.R")
source("../../functions/import_data.R")

d <- import_data("clarke2022qjep")

m13te <- readRDS("../1_test_old_model_old_data/scratch/clarke2022qjep_test_1_3.model")
m13tr <- readRDS("../1_test_old_model_old_data/scratch/clarke2022qjep_train_1_3.model")

pred13 <- summarise_postpred(list(training = m13tr, testing = m13te), 
                             d,  multi_level = TRUE, get_sim = FALSE)

rm(m13te, m13tr)

m12te <- readRDS("../1_test_old_model_old_data/scratch/clarke2022qjep_test_1_2.model")
m12tr <- readRDS("../1_test_old_model_old_data/scratch/clarke2022qjep_train_1_2.model")

pred12 <- summarise_postpred(list(training = m12tr, testing = m12te), 
                             d,  multi_level = TRUE, get_sim = FALSE)

rm(m12te, m12tr)

acc12 <- compute_acc(pred12$acc)
acc13 <- compute_acc(pred13$acc)

rm(pred10, pred13)

bind_rows(acc12 %>% mutate(model = "1.2"),
          acc13 %>% mutate(model = "1.3")) %>%
  ggplot(aes(x = found, y = accuracy, fill = model,
                  ymin = .lower, ymax = .upper, group = interaction(model, .width))) +
  geom_lineribbon(alpha = 0.75) + 
  facet_grid(condition~split)
