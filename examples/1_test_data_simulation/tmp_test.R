library(tidyverse)


options(mc.cores = 8, 
        digits = 2)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

source("../../functions/sim_foraging_data.R")
source("../../functions/plot_data.R")


trl = 1 
n_item_class = 4
n_item_per_class = 20
item_labels = c("A", "B", "d1", "d2")


## Proximity Bias Only


item_class_weights = c(0.5, 0.5, 0, 0)

b_stick  <- 0  # no stick or switching
b_memory <- 0 # no memory of previous weights

rho_delta <- 5# proximity tuning
rho_psi   <- 0  # rel. direction tuning
abs_dir_tuning <- list(kappa = rep(10, 4), theta = c(3, 3, 3, 3))

# initial bias params
inital_sel_params <- tibble(
  a1x = 2,
  b1x = 2,
  a2x = 1,
  b2x = 10,
  a1y = 2,
  b1y = 2,
  a2y = 10,
  b2y = 1) 

init_sel_lambda <- 0

items <- tibble(x = seq(0, 1, 0.2), y = x) %>% modelr::data_grid(x,y ) %>%
  mutate(item_class = 1, id = 1:n())

d <- sim_foraging_trial(trl = 1, 
                        n_item_class, n_item_per_class, 
                        item_class_weights,
                        item_labels,
                        b_stick, 
                        rho_delta, 
                        rho_psi,
                        abs_dir_tuning,
                        b_memory,
                        inital_sel_params = inital_sel_params,
                        init_sel_lambda = init_sel_lambda,
                        items = items,
                        dev_output = TRUE)

# plot absolute tuning curve
ddir <- tibble(phi = seq(1, 360),
            z = compute_all_von_mises(phi=phi, abs_dir_tuning$theta, abs_dir_tuning$kappa))

ggplot(ddir, aes(phi, z)) + geom_path()


plot_a_trial(d$stim, d$found, "delta")

