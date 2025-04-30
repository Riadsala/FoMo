# minimum working example for FoMo 1.3 and higher (includes direction weight plots etc)

# 1) simulate some data (one person, one condition, multiple trials)
# 2) prep data
# 3) fit model
# 4) extract, summarise and plot posterior
# 5) (not yet implemented for single-level models) generated quantities (see multi-level example)

library(tidyverse)
library(cmdstanr)


source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")


options(mc.cores = 4)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

######################################################################
# lets simulate some data - this is for the TEST dataset
######################################################################
n_trials_per_cond <- 10
n_items_per_class <- 14

stimuli_params <- list(
  n_item_class = 4,
  n_item_per_class = c(n_items_per_class, n_items_per_class, n_items_per_class, n_items_per_class), 
  item_labels = c("a", "b", "d1", "d2"))

foraging_params <- list(
  b_a = 0.5, 
  b_s = 1, 
  rho_delta = 1,
  rho_psi = 0)

absdir_params = list(
  kappa = rep(25, 4), theta = c(5, 3, 5, 3))

# # plot absolute tuning curve
# d <- tibble(phi = seq(0, 2*pi, 0.01),
#             z = compute_all_von_mises(phi=phi, absdir_params$theta, absdir_params$kappa))
# 
# ggplot(d, aes(phi, z)) + geom_path()


# simulate data
d <- sim_foraging_multiple_trials(person = 1, 
                                  n_trials_per_cond = n_trials_per_cond,
                                  condition = "test",
                                  sp = stimuli_params,
                                  fp = foraging_params,
                                  adp = absdir_params, 
                                  isp = "off",
                                  dev_output = TRUE)

# check empirical distribution
plot_a_trial(d$stim, d$found, 1, segLabel = "phi")

######################################################################
# prep data
######################################################################

modelver <- "1.3"
modelver_str <- str_replace(modelver, "\\.", "_" )

dl <- prep_data_for_stan(d)
dl <- add_priors_to_d_list(dl, modelver = modelver)

######################################################################
# fit model
######################################################################

iter = 500
mod <- cmdstan_model(paste0("../../models/simple/FoMo", modelver_str, ".stan"))
m <- mod$sample(data = dl, 
                iter_warmup  = iter, iter_sampling = iter)

######################################################################
# check posterior
######################################################################

m$summary() %>% knitr::kable(digits = 2) 
post <- extract_post(m, d)

# plot core parameters
plot_model_fixed(post, foraging_params)

# plot directional component
plot_model_theta(post)

post$theta %>% ggplot(aes(comp, log(theta))) +
  stat_interval()
######################################################################
# diagnostics
######################################################################

bayesplot::mcmc_trace(m$draws(), pars = "lp__")
