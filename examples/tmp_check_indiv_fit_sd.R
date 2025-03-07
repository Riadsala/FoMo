library(tidyverse)
library(cmdstanr)
library(patchwork)
library(tidybayes)

source("../functions/import_data.R")
source("../functions/prep_data.R")
source("../functions/compute_summary_stats.R")
source("../functions/plot_model.R")
source("../functions/plot_data.R")
source("../functions/post_functions.R")
source("../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

model_ver <- "1_0"
dataset <- "hughes2024rsos"

# read in data
d <- import_data(dataset)

sf <- "1_fit_models/scratch"
folder <- paste0(sf, "/post/", dataset, "/")



m <- readRDS(paste0("1_fit_models/scratch/models/", dataset, "/train", model_ver, ".model"))
post <- extract_post(m, d)
plot_model_random(post)

# plot variances
prior <- tibble(value = rexp(1000, 1))

# calculate sd from the post indiv diff estimates
post$random %>%
  pivot_longer(c(u_a, u_stick, u_delta), names_to = "param") %>%
  group_by(.draw, param, condition) %>%
  summarise(sd = sd(value)) -> calc_from_u

post$variances %>% 
  mutate(param = str_replace_all(param, "b|rho", "u")) %>%
  ggplot(aes(value)) +
  geom_density(data = prior, fill = "grey") +
  geom_density(aes(fill = condition), alpha = 0.5) + 
  stat_intervalh(data = calc_from_u, aes(colour = condition, x = sd), y = 0) + 
  facet_wrap(~param)


# now repeat for the theta
prior <- tibble(value = rexp(1000, 1))

post$utheta %>%
  group_by(.draw, comp, condition) %>%
  summarise(sd = sd(log(theta))) -> calc_from_u

post$theta %>%
  ggplot(aes(sigma)) +
  geom_density(data = prior, aes(value), fill = "grey") +
 # geom_density(aes(fill = condition), alpha = 0.5) + 
  stat_intervalh(data = calc_from_u, aes(colour = condition, x = sd), y = 0) + 
  facet_wrap(~comp)

