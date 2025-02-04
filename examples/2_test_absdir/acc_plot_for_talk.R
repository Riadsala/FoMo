library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(patchwork)

source("../../functions/post_functions.R")
source("../../functions/import_data.R")
source("../../functions/plot_data.R")

d <- import_data("clarke2022qjep")

m13te <- readRDS("scratch/clarke2022qjep_test_1_3.model")
m13tr <- readRDS("scratch/clarke2022qjep_train_1_3.model")

pred13 <- summarise_postpred(list(training = m13tr, testing = m13te), 
                             d,  multi_level = TRUE, get_sim = FALSE)

rm(m13te, m13tr)

m12te <- readRDS("scratch/clarke2022qjep_test_1_2.model")
m12tr <- readRDS("scratch/clarke2022qjep_train_1_2.model")

pred12 <- summarise_postpred(list(training = m12tr, testing = m12te), 
                             d,  multi_level = TRUE, get_sim = FALSE)

rm(m12te, m12tr)

acc12 <- compute_acc(pred12$acc)
acc13 <- compute_acc(pred13$acc)

rm(pred10, pred13)


chance  = tibble(found = 1:40,
                 accuracy = 1/(41-found))

bind_rows(acc12 %>% mutate(model = "1.2"),
          acc13 %>% mutate(model = "1.3")) %>%
  ggplot(aes(x = found, y = accuracy)) +
  geom_lineribbon(alpha = 0.75, aes(fill = model, colour = model,
                                    ymin = .lower, ymax = .upper, group = interaction(model, .width))) + 
  facet_grid(condition~split) + 
  geom_path(data = chance, aes(found, accuracy), linetype = 2) +
  theme_bw()

ggsave("acc.png", width = 4,  height = 3)

acc13 <- compute_acc(pred13$acc, compute_hpdi = FALSE)
acc12 <- compute_acc(pred12$acc, compute_hpdi = FALSE)


bind_rows(acc12 %>% mutate(model = "1.2"),
          acc13 %>% mutate(model = "1.3")) %>%
  group_by(split, model, condition) %>%
  filter(found < 30) %>%
  median_hdci(accuracy)


post <- extract_post(m13tr, d)

post$utheta %>% 
  group_by(person, comp, condition) %>%
  summarise(phi = unique(phi),
            theta = median(theta)) %>%
  ggplot(aes(phi, theta, fill = condition)) + 
  geom_col(position = position_identity(), alpha = 0.5)  +
  facet_wrap(~person, nrow = 6) +
  theme_bw() +
  ggthemes::scale_fill_ptol()

ggsave("indiv_dir_plots.png", width = 8, height = 8)


set.seed(2025)
plot_a_trial(d$stim, d$found, sample(filter(d$found, person == 2)$trial, 1)) +
  ggtitle(paste("person", 2)) -> plt10

plot_a_trial(d$stim, d$found, sample(filter(d$found, person == 34)$trial, 1)) +
  ggtitle(paste("person", 34)) -> plt20

plot_a_trial(d$stim, d$found, sample(filter(d$found, person == 6)$trial, 1)) +
  ggtitle(paste("person", 6)) -> plt40

plt10 / plt20 / plt40

ggsave("examples.png", width = 6, height = 12)

post$theta %>%s
  ggplot(aes(phi, theta, fill = condition, colour = condition)) + 
  stat_slabinterval(alpha = 0.5)+
  theme_bw() +
  ggthemes::scale_fill_ptol()+
  ggthemes::scale_colour_ptol() +
  scale_x_discrete("direction") + 
  theme(legend.position = "bottom")

ggsave("group_avg.png", width = 6, height = 8)
