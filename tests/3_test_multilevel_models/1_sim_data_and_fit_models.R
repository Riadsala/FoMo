library(tidyverse)
library(cmdstanr)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")

options(mc.cores = 8)

######################################################################
# lets simulate some data 
######################################################################

experiment_params <- list(n_people = 4, 
                          n_conditions = 2,
                          condition_labels = c("A", "B"),
                          n_trials_per_cond = 5)

stimuli_params <- list(n_item_class = 4,
                       n_item_per_class = c(10, 10, 10, 10), 
                       item_labels = c("a", "b", "d1", "d2"))

foraging_params <- list(b_a = c(0, 0.25), 
                        b_s = c(0.5, 2.5), 
                        rho_delta = c(1, 0.8),
                        rho_psi = c(-0.5, 0.5))

variance_params <- list(b_a = c(0.1, 0.1), 
                        b_s = c(0.1, 0.1), 
                        rho_delta = c(0.1, 0.1),
                        rho_psi = c(0.1, 0.1))

absdir_params <- list(
  kappa = rep(20, 4), theta = c(5, 1, 5, 1))

initsel_params <- "off"

params <- list(e = experiment_params,
               s = stimuli_params,
               f = foraging_params,
               v = variance_params,
               a = absdir_params,
               i = initsel_params)

d <- sim_foraging_people(params, filename = "data_for_qmd") 


######################################################################
# fit model, generated quantities
######################################################################

modelver <- c("1.3")
modelname <- "qmd_model"

for (mv in modelver) { 
  
  modelver_str <- str_replace(mv, "\\.", "_" )
  
  dl <- prep_data_for_stan(d)
  dl <- add_priors_to_d_list(dl, modelver = mv)

  iter = 500
  mod <- cmdstan_model(paste0("../../models/multi_level/FoMo", modelver_str, ".stan"),
                      cpp_options = list(stan_threads = TRUE))

  m <- mod$sample(data = dl, 
                  iter_warmup  = iter, iter_sampling = iter,
                  chains = 4, 
                  parallel_chains = 8,
                  threads_per_chain = 2)

  m$save_object(paste0("scratch/fit/", modelname, "_", modelver_str, ".model"))
  
  # generated quantities
  
  draws_matrix <- posterior::as_draws_matrix(m$draws())
  idx <- sample(nrow(draws_matrix), 100) #iter_genquant
  
  mod_sim <- cmdstan_model(paste0("../../models/simulate/FoMo", modelver_str, ".stan"))
  
  p <- mod_sim$generate_quantities(fitted_params = draws_matrix[idx,],
                                   data = dl,
                                   seed = 123,
                                   output_dir = "scratch/sim",
                                   output_basename = paste(modelname, "_", modelver_str, sep=""))
  
  p$save_object(paste0("scratch/sim/", modelname, "_", modelver_str, ".model"))
  

}

######################################################################
# extract predictions
######################################################################

# compute empirical run statistics
rl <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs = mean(n_runs),
            mean_bestr = mean(best_r),
            mean_pao = mean(pao),
            .groups = "drop") %>%
  mutate(z = "observed")

# compute empirical run statistics
iisv <- get_iisv_over_trials(d$found) %>%
  mutate(z = "observed")
