
library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")

options(mc.cores = 4)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

d <- import_data('tagu2022cog')

pp <- 4
  
d_one_person <- filter_one_person(d, pp) 

dloo <- readRDS("scratch/person4loo11.rds")

d_one_person$found %>% 
  mutate(high_pareto = if_else(
  dloo$diagnostics$pareto_k > 0.7,
  "over 0.7", "below 0.7")) -> d_pareto

d_pareto %>% group_by(trial) %>%
  summarise(p = sum(high_pareto == "over 0.7"))

d_pareto %>% 
  filter(trial == 10) %>% 
ggplot(aes(x, y, colour = high_pareto)) + 
  
  geom_path()
