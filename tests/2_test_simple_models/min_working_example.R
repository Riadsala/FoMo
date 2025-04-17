# minimum working example
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
n_trials_per_cond <- 50

stimuli_params <- list(n_item_class = 4,
                       n_item_per_class = c(20, 20, 10, 10), 
                       item_labels = c("a", "b", "d1", "d2"))

foraging_params <- list(
  b_a = 0.5, 
  b_s = 2, 
  rho_delta = 1.5,
  rho_psi = -0.5)

d <- sim_foraging_multiple_trials(person = 1, 
                                  n_trials_per_cond = n_trials_per_cond,
                                  condition = "test",
                                  sp = stimuli_params,
                                  fp = foraging_params,
                                  adp = "off", isp = "off")

######################################################################
# prep data
######################################################################

modelver <- "1.0"
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
plot_model_fixed(post, foraging_params)

######################################################################
# diagnostics
######################################################################

bayesplot::mcmc_trace(m$draws(), pars = "lp__")
