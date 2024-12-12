library(tidyverse)
library(cmdstanr)
library(tidybayes)

source("../../functions/post_functions.R")
source("../../functions/import_data.R")

d <- import_data("hughes2024rsos")

m <- readRDS("scratch/fomo13_hughes2024.model")

post <- extract_post(m, d)


post$utheta %>% ggplot(aes(factor(phi), theta), colour = factor(person)) +
  stat_interval(alpha = 0.2, 
                position = position_dodge())


compute_von_mises <- function(x, .draw, person, condition, phi, theta, kappa, comp) {
  
  z <- theta * exp(kappa * cos(phi-x)) / (2*pi*besselI(kappa,0))
  
  dout <- tibble(.draw = .draw,
                 person = person,
                 condition = condition,
                 x = x, 
                 z = z,
                 comp = comp)
  return(dout)
  
}

post$utheta %>% mutate(kappa = 10) %>%
  filter(.draw < 500) %>%
  select(.draw, person, condition, phi, theta, kappa, comp) %>%
  pmap_df(compute_von_mises, x = seq(0, 2*pi, 0.1)) %>%
  group_by(.draw, person, condition, x) %>%
  summarise(weight = sum(z) + 1) -> post_weights

post_weights %>%
  filter(.draw < 50) %>%
  ggplot(aes(x, weight, 
             color = condition, 
             group = interaction(condition, .draw))) +
  geom_path(alpha = 0.1) + 
  facet_wrap(~person, scales = "free") +
  coord_polar() +
  theme_minimal()

trl <- 2

pp <- 10

ggplot() +
  geom_point(data = d$stim %>% 
               filter(person == pp, trial_p == trl,), 
             aes(x, y, colour = factor(item_class))) +
  geom_path(data = d$found %>% 
               filter(person == pp, trial_p == trl), 
             aes(x, y)) + 
  facet_wrap(~condition)


pred13 <- summarise_postpred(m, d,  multi_level = FALSE, get_sim = FALSE)


m12 <- readRDS("../1_test_old_model_old_data/scratch/hughes2024rsos_test_1_2.model")
pred12 <- summarise_postpred(m12, d,  multi_level = FALSE, get_sim = FALSE)


pred13$acc %>%
  group_by(found, condition, .draw) %>% 
  summarise(acc = mean(model_correct)) %>%
  median_hdci(acc) -> acc95


ggplot(acc95, aes(found, 
                  ymin = .lower, y = acc, ymax = .upper, fill = condition)) + 
  geom_lineribbon(alpha = 0.5) +
  coord_cartesian(ylim = c(0.4, 0.55))

