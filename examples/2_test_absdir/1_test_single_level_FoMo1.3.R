library(tidyverse)
library(patchwork)
library(cmdstanr)

### Initial simple example
# Test model on single level data

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")
source("../../functions/prep_data.R")

d <- import_data("hughes2024rsos")

person <- 2

d$stim  %>% filter(person == {{person}}) -> d$stim
d$found %>% filter(person == {{person}}) -> d$found

d$stim <- fix_person_and_trial(d$stim)
d$found <- fix_person_and_trial(d$found)

d_list <- prep_data_for_stan(d$found, d$stim, c("spatial", "item_class"))

d_list$prior_mu_b_a <- 0
d_list$prior_sd_b_a <- 0.5
d_list$prior_mu_b_stick <- 0
d_list$prior_sd_b_stick <- 1
d_list$prior_mu_rho_delta <- 15
d_list$prior_sd_rho_delta <- 5
d_list$prior_mu_rho_psi <- 0
d_list$prior_sd_rho_psi <- 1
d_list$prior_theta_lambda <- 0.1
d_list$n_trials_to_sim <- 3

d_list$kappa <- 10

iter = 500


run_FoMo <- function(sl, model_ver) { 
  
  model_fn <- paste0("../../models/simple/FoMo", model_ver, ".stan")

  mod <- cmdstan_model(model_fn, 
                       cpp_options = list(stan_threads = TRUE))
  
  # run model
  m <- mod$sample(data = d_list, 
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 100, 
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3)
  
  return(m)
}

# run model
m13 <- run_FoMo(d_list, "1_3")
m12 <- run_FoMo(d_list, "1_2")

# extract post
post <- extract_post(m, d, multi_level = FALSE, absdir = TRUE)

# plot model
plot_model_fixed(post)

ggplot(post$absdir, aes(theta, fill = factor(comp))) + geom_density(alpha = 0.4) 


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
pred12 <- summarise_postpred(m12, d,  multi_level = FALSE)
pred13 <- summarise_postpred(m13, d,  multi_level = FALSE)

bind_rows(pred12$acc %>% mutate(model = "1.2"),
          pred13$acc %>% mutate(model = "1.3")) -> pred


n_targets <- max((pred$found))

baseline <- tibble(found = 1:n_targets, accuracy = 1/((n_targets + 1) - found))

pred %>% 
  group_by(model, condition, found, .draw) %>%
  summarise(accuracy = mean(model_correct, .groups ="last")) %>%
  median_hdci(accuracy) %>%
  ggplot(aes(found, accuracy)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model), alpha = 0.5) + 
  geom_path() + 
  geom_path(data = baseline, linetype = 2) +
  facet_wrap(~condition)


# plot comparison between a real and simulated trial

pltreal <- plot_a_trial(d$stim, d$found, 1)
pltsim <- plot_a_trial(d$stim, pred$sim %>% filter(.draw == 1), trial = 1)

pltreal + pltsim

