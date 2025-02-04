library(tidyverse)
library(cmdstanr)
library(tidybayes)
library(patchwork)

source("../../functions/post_functions.R")
source("../../functions/import_data.R")
source("../../functions/plot_data.R")

d <- import_data("clarke2022qjep")


mte <- readRDS("scratch/clarke2022qjep_test_1_0.model")
mtr <- readRDS("scratch/clarke2022qjep_train_1_0.model")

pred <- summarise_postpred(list(training = mtr, testing = mte), 
                             d,  multi_level = TRUE, get_sim = FALSE)



acc <- compute_acc(pred$acc)

chance  = tibble(found = 1:40,
                 accuracy = 1/(41-found))

acc %>%
  ggplot(aes(x = found, y = accuracy)) +
  geom_lineribbon(alpha = 0.75, aes(fill = condition, colour = condition,
                                    ymin = .lower, ymax = .upper, group = interaction(condition, .width))) + 
  facet_grid(.~split) + 
  geom_path(data = chance, aes(found, accuracy), linetype = 2) +
  theme_bw() + 
  ggthemes::scale_fill_ptol() + 
  ggthemes::scale_colour_ptol()

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
