# ignoring multilevels
# does everything work?

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
n_trials_per_cond <- 25

stimuli_params <- list(n_item_class = 4,
                       n_item_per_class = c(10, 10, 10, 10), 
                       item_labels = c("a", "b", "d1", "d2"))

foraging_params <- list(
  b_a = 0, 
  b_s = 0, 
  rho_delta = 1,
  rho_psi = 0)

d <- sim_foraging_multiple_trials(person = 1, 
                                  n_trials_per_cond = n_trials_per_cond,
                                  condition = "test",
                                  sp = stimuli_params,
                                  fp = foraging_params,
                                  adp = "off", isp = "off")

######################################################################
# lets simulate some data - this is for the TEST MULTICOND dataset
######################################################################

experiment_params <- list(n_people = 1, 
                          n_conditions = 2,
                          condition_labels = c("A", "B"),
                          n_trials_per_cond = 25)

stimuli_params <- list(n_item_class = 4,
                       n_item_per_class = c(10, 10, 10, 10), 
                       item_labels = c("a", "b", "d1", "d2"))

foraging_params <- list(b_a = c(0, 1), 
                        b_s = c(2, 0), 
                        rho_delta = c(1, 1),
                        rho_psi = c(-1, 1))

variance_params <- list(b_a = c(0.1, 0.1), 
                        b_s = c(0.1, 0.1), 
                        rho_delta = c(0.1, 0.1),
                        rho_psi = c(0.1, 0.1))

absdir_params <- "off"
initsel_params <- "off"

params <- list(e = experiment_params,
               s = stimuli_params,
               f = foraging_params,
               v = variance_params,
               a = absdir_params,
               i = initsel_params)


d <- sim_foraging_people(params) 



######################################################################
# prep data
######################################################################

modelver <- "1.3"
modelname <- "test_multicond"
modelver_str <- str_replace(modelver, "\\.", "_" )

saveRDS(d, paste0("scratch/", modelname, "_data.RDS"))


dl <- prep_data_for_stan(d)
dl <- add_priors_to_d_list(dl, modelver = modelver)

######################################################################
# fit model
######################################################################

iter = 500
mod <- cmdstan_model(paste0("../../models/simple/FoMo", modelver_str, ".stan"))
m <- mod$sample(data = dl, 
                iter_warmup  = iter, iter_sampling = iter)

m$save_object(paste0("scratch/", modelname, "_", modelver_str, ".model"))

######################################################################
# check posterior
######################################################################

m$summary() %>% knitr::kable(digits = 2) 
post <- extract_post(m, d)
plot_model_fixed(post, foraging_params)
