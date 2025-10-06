library(tidyverse)
library(cmdstanr)
library(patchwork)
library(tidybayes)

source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/post_functions.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

dataset <- "bhat2025"

# read in data
d <- import_data(dataset)

sf <- "../1_fit_models/scratch"
folder <- paste0(sf, "/post/", dataset, "/")

# compare accuract across model versions

plt_acc_2 <- plot_models_accuracy(dataset, sf)

plt_acc_2


# comparing 1.3 with 1.0

v1 <- "1_0_k20"
v2 <- "1_3_k20"

plt_acc_3 <- plot_model_accuracy_comparison(dataset, v1, v2, sf) 

plt_acc_3

# fixed effect plots

model_ver <- "1_3"


m <- readRDS(paste0("1_fit_models/scratch/models/", dataset, "/fit/", model_ver, ".model"))

post <- extract_post(m, d)
post_plt <- plot_model_fixed(post)
post_plt

# need to do a little rejigging because of grid offsets
post$theta <- post$theta %>%
  mutate(phi = ifelse(str_detect(condition, "rot"), phi + pi/4, phi))

post_theta <- plot_model_theta(post, nrow = 1)
post_theta

# a table?

datasets <- "bhat2025"#c("kristjansson2014plos", "tagu2022cog", "hughes2024rsos", "clarke2022qjep")

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
  group_by(dataset, model_ver, condition, .draw) %>%
  summarise(accuracy = mean(model_correct), .groups = "drop_last") %>%
  summarise(accuracy = mean(accuracy), .groups = "drop") %>%
  mutate(model_ver = str_replace(model_ver, "_", ".")) %>%
  pivot_wider(names_from = model_ver, values_from = accuracy) %>%
  full_join(d_chance, by = join_by(dataset)) %>%
  knitr::kable()

# plot a trial

ds <- d$stim
df <- d$found

good_grid <- plot_a_trial(ds, df, 2) # person 1, conj grid, for whom 1.3 improves things a lot (accuracy goes from quite low to quite middling)
plot_a_trial(ds, df, 15) # person 1, feat grid, 1.3 improves things a lot (accuracy goes from lowish to middling)
plot_a_trial(ds, df, 50) # person 3, conj grid, for whom 1.3 improves things (accuracy stays relatively low)

plot_a_trial(ds, df, 43) # person 2, feat rot, 1.3 actually makes things slightly worse but accuracy is v good
plot_a_trial(ds, df, 133) # person 6, feat grid, high accuracy, goes up a little in 1.3 but I think mostly it's just proximity
plot_a_trial(ds, df, 139) # person 6, feat rot, high accuracy - again, I think mostly proximity (but does go up a little with 1.3)

plot_a_trial(ds, df, 265) # person 12, conj grid, who doesn't show any difference really with 1.3 (quite low accuracy overall, mostly because I don't think they're very proximity based)
plot_a_trial(ds, df, 638) # person 27, feat grid v low (I think because they are doing things in a conjunction-y way? Or changing a lot between trials?)
plot_a_trial(ds, df, 549) # person 23, feat rot, quite high accuracy but goes down a little for 1.3 (does seem to be going up and down a bit)

bad_grid <- plot_a_trial(ds, df, 647) # person 28, conj grid, not doing a very good job (and 1.3 not helping)
plot_a_trial(ds, df, 653) # person 28, conj rot, not that great and 1.3 not helping. Search does look a bit disorganised and not particularly on any grid.
plot_a_trial(ds, df, 659) # person 28, feat grid, doing  better (1.3 helps a bit?)
plot_a_trial(ds, df, 665) # person 28, feat rot, doing better (1.3 helps a bit?)

good_grid + bad_grid

ggsave('example_trials.png', bg = "transparent")
