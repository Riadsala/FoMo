library(tidyverse)
library(cmdstanr)
library(loo)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")

# create scratch folder if it does not yet exist
if (!file.exists("scratch"))  dir.create("scratch")

item_class_weights = list(c(0.7, 0.3, 0, 0))

b_stick = 1

rho_delta = 20 #1 - set to 1 if doing relative proximity = TRUE, 20 if relative proximity = FALSE
sd_rho_delta = 1

rho_psi = 2

abs_dir_tuning = list(kappa = rep(10, 4), theta = rep(1, 4))

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

d <- sim_foraging_people(n_people = 8,
                         n_conditions = 1,
                         cond_lab = c("simple test"),
                         n_trials_per_cond = 10,
                         n_item_class = 2, n_item_per_class = 20,
                         item_class_weights, sd_bA = 0.2,
                         b_stick = b_stick, sd_b_stick = 1,
                         rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                         rho_psi = rho_psi, sd_rho_psi = 0.25,
                         abs_dir_tuning = abs_dir_tuning,
                         inital_sel_params = inital_sel_params,
                         rel_proximity = FALSE) 


d_list <- prep_data_for_stan(d$found, d$stim, model_components = c("spatial", "item_class"), n_trials_to_sim = 2)
d_list <- add_priors_to_d_list(d_list, modelver = "1.3")

mod <- cmdstan_model(paste0("../../models/multi_level/FoMo1_3.stan"))

m <- mod$sample(data = d_list, 
                chains = 4, parallel_chains = 4, threads = 4,
                refresh = 10, 
                iter_warmup = 500, iter_sampling = 500,
                sig_figs = 3)

############################################################
############################################################

# extract post
post <- extract_post(m, d, multi_level = FALSE)

# plot model
plot_model_fixed(post, gt = list(b_a = plogis(item_class_weights[[1]]),
                                 b_stick = b_stick,
                                 rho_delta = rho_delta,
                                 rho_psi = rho_psi))

ggplot(post$absdir, aes(theta, fill = factor(comp))) + geom_density(alpha = 0.4) +
  geom_vline(data = tibble(comp = 1:4, x = abs_dir_tuning$theta), aes(xintercept=x, colour = factor(comp)))

# plot von mises distribution
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


# check predictions and accuracy
pred <- summarise_postpred(m, d, draw_sample_frac = 0.1, multi_level = FALSE, get_sim = TRUE)

plot_model_accuracy(pred)


