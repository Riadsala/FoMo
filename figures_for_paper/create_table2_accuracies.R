library(tidyverse)
library(cmdstanr)

source("../functions/import_data.R")
source("../functions/prep_data.R")
source("../functions/compute_summary_stats.R")
source("../functions/plot_model.R")
source("../functions/plot_data.R")
source("../functions/post_functions.R")
source("../functions/sim_foraging_data.R")

sf <- "../examples/1_fit_models/scratch"
datasets <- c("kristjansson2014plos", "tagu2022cog", "hughes2024rsos", "clarke2022qjep")

d_acc <- tibble()
d_chance <- tibble()

for (ds in datasets) {
  # find list of models
  files <- dir(paste0(sf, "/post/", ds))
  files <- files[str_detect(files,"pred")]
  
  d <- import_data(ds)
  n_targets <- max((d$found$found))
  baseline <- tibble(found = 1:n_targets, accuracy = 1/((n_targets + 1) - found))
  d_chance <- bind_rows(d_chance, tibble(dataset = ds, chance = mean(baseline$accuracy)))
  rm(d, baseline, n_targets)
  
  
  for (ff in files) {
    
    a <- readRDS(paste0(sf, "/post/", ds, "/", ff))
    a$itemwise$model_ver <- a[[5]]
    d_acc <- bind_rows(d_acc, a$itemwise %>% mutate(dataset = ds))
    
  }
}

# create our table
d_acc %>% 
  group_by(dataset, model_ver, .draw) %>%
  summarise(accuracy = mean(model_correct), .groups = "drop_last") %>%
  summarise(accuracy = mean(accuracy), .groups = "drop") %>%
  mutate(model_ver = str_replace(model_ver, "_", ".")) %>%
  pivot_wider(names_from = model_ver, values_from = accuracy) %>%
  full_join(d_chance, by = join_by(dataset)) %>%
  knitr::kable()
