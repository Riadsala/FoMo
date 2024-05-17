library(tidyverse)
library(patchwork)
library(cmdstanr)

### Initial simple example
# Test model on single level data

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

source("../../functions/sim_foraging_data.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")
source("../../functions/prep_data.R")

n_trials_per_cond <- 50

n_item_class <- 2
n_item_per_class <- 20
item_class_weights = c(0.7, 0.3, 0, 0)
b_stick = 2
b_memory = 0

abs_dir_tuning = list(kappa = rep(10, 4), theta = c(0, 0, 0, 0))
rho_delta = 10
rho_psi = 5

d <- sim_foraging_multiple_trials(person = 1, 
                                  condition = "test",
                                  n_item_class =  n_item_class, n_item_per_class = n_item_per_class,
                                  item_class_weights = item_class_weights, item_labels = item_labels,
                                  b_stick = b_stick, 
                                  rho_delta = rho_delta, 
                                  rho_psi = rho_psi, 
                                  abs_dir_tuning = abs_dir_tuning,
                                  b_memory = b_memory,
                                  inital_sel_params = inital_sel_params,
                                  init_sel_lambda = init_sel_lambda)

iter = 500
mod <- cmdstan_model("../../models/simple/FoMo1_2.stan", 
                     cpp_options = list(stan_threads = TRUE))

d_list <- prep_data_for_stan(d$found, d$stim, c("spatial", "item_class"))

# add priors to list
d_list$prior_mu_b_a <- 0
d_list$prior_sd_b_a <- 0.5
d_list$prior_mu_b_stick <- 0
d_list$prior_sd_b_stick <- 1
d_list$prior_mu_rho_delta <- 15
d_list$prior_sd_rho_delta <- 5
d_list$prior_mu_rho_psi <- 0
d_list$prior_sd_rho_psi <- 1
d_list$n_trials_to_sim <- 10

d_list$kappa <- 10

# run model
m <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 100, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

# extract post
post <- extract_post(m, d, multi_level = FALSE, absdir = TRUE)

# plot model
plot_model_fixed(post, gt = list(b_a = plogis(item_class_weights[[1]]),
                                 b_stick = b_stick,
                                 rho_delta = rho_delta,
                                 rho_psi = rho_psi))

ggplot(post$absdir, aes(theta, fill = factor(comp))) + geom_density(alpha = 0.4) +
  geom_vline(data = tibble(comp = 1:4, x = abs_dir_tuning$theta), aes(xintercept=x, colour = factor(comp)))


compute_von_mises <- function(x, .draw, phi, theta, kappa, comp) {
  
  z <- theta * exp(kappa * cos(phi-x)) / (2*pi*besselI(kappa,0))
  
  dout <- tibble(.draw = .draw,
                 x = x, 
                 z = z,
                 comp = comp)
  return(dout)
  
}

post$absdir %>% mutate(kappa = 10) %>%
  select(.draw, phi, theta, kappa, comp) %>%
  pmap_df(compute_von_mises, x = seq(0, 2*pi, 0.01)) %>%
  group_by(.draw, x) %>%
  summarise(weight = sum(z) + 1) -> post_weights

post_weights %>% 
  group_by(x) %>%
  median_hdci(weight, .width = c(0.53, 0.97)) %>%
  ungroup() %>%
  ggplot(aes(x, weight, ymin = .lower, ymax = .upper, group = .width)) + 
  geom_ribbon(alpha = 0.5, fill = "purple")


# check predictions
pred <- summarise_postpred(m, d)

plot_model_accuracy(pred)

# plot comparison between a real and simulated trial

pltreal <- plot_a_trial(d$stim, d$found, 1)
pltsim <- plot_a_trial(d$stim, pred$sim %>% filter(.draw == 1), trial = 1)

pltreal + pltsim

####### Now check FoMo1.1 works

mod <- cmdstan_model("../../models/simple/FoMo1_1.stan", 
                     cpp_options = list(stan_threads = TRUE))
# run model
m <- mod$sample(data = d_list, 
                chains = 4, parallel_chains = 4, threads = 4,
                refresh = 100, 
                iter_warmup = iter, iter_sampling = iter,
                sig_figs = 3)
