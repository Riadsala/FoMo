library(tidyverse)
library(patchwork)
library(circular)


options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

library(tidyverse)
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
theme_set(theme_bw())

theta <- c(15, 4, 8, 6)
mu <- c(0, pi/2, pi, 3*pi/2)
kappa <- 10


phi <- seq(-pi, pi, 0.01) 

comp_vm <- function(k, theta, mu, kappa, phi) {
  
  tibble(k = k, 
         phi = phi,
         w = theta * dvonmises(phi, mu, kappa))
  
}

d <- tibble(k = 1:4, theta = theta, mu = mu, kappa = kappa) %>%
  pmap_df(comp_vm,  phi) %>%
  mutate(k = factor(k))

ggplot(d, aes(phi, w, colour = k)) + 
  geom_path(linewidth = 2) +
  geom_hline(yintercept = 0, linewidth = 2) + 
  scale_x_continuous(limits = c(-pi, pi), expand = c(0,0)) +
  theme(legend.position = "none") -> plt1

d %>% group_by(phi) %>%
  summarise(w = sum(w)) %>%
  ggplot(aes(phi, w)) + 
  geom_path(linewidth = 2) +
  scale_x_continuous(limits = c(-pi, pi), expand = c(0,0)) -> plt2

plt1 / plt2
  
ggsave("figs/figure4_vm.pdf", width = 6, height = 4)

