
library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/plot_model.R")
source("../../functions/plot_data.R")
source("../../functions/post_functions.R")
source("../../functions/sim_foraging_data.R")

options(mc.cores = 4)

# set global ggplot theme
theme_set(ggthemes::theme_tufte())

d <- import_data('tagu2022cog')

mod11 <- cmdstan_model("../../models/simple/FoMo1_1.stan", 
                     cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)

mod12 <- cmdstan_model("../../models/simple/FoMo1_2.stan", 
                     cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)

iter = 500


for (pp in 1:24) {

  d_one_person <- filter_one_person(d, pp) 
  
  d_list <- prep_data_for_stan(d_one_person$found, d_one_person$stim, c("spatial", "item_class"))
  d_list <- add_priors_to_d_list(d_list)
  d_list$n_trials_to_sim <- 1
  
  filename <- paste0("scratch/person", pp)

  # model 1.1
  fit <- mod11$sample(data = d_list, 
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 100, 
                    init = 1,
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3)
  
 fit$save_object(paste0(filename, "_11.rds"))
  
  # model 1.2
  d_list$prior_theta_lambda <- 10
  d_list$kappa <- 10
  
  fit <- mod12$sample(data = d_list, 
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 100, 
                    init = 1,
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3)
  
  fit$save_object(paste0(filename, "_12.rds"))
  
}

for (pp in 1:24) {
  
  filename <- paste0("scratch/person", pp)
  
  m11 <- readRDS(paste0(filename, "_11.rds"))
  m12 <- readRDS(paste0(filename, "_12.rds"))
  
  loo11 <- m11$loo()
  loo12 <- m12$loo()
  
  saveRDS(loo11, paste0(filename, "loo11.rds"))
  saveRDS(loo12, paste0(filename, "loo12.rds"))
  
}

w <- tibble()

for (pp in 1:24) {
  
  filename <- paste0("scratch/person", pp)
  
  loo11 <- readRDS(paste0(filename, "loo11.rds"))
  loo12 <- readRDS(paste0(filename, "loo12.rds"))
  
  t <- loo::loo_model_weights(list("FoMo 1.1" = loo11, "FoMo 1.2" = loo12))
  
  w <- bind_rows(w, tibble(
    person = pp, FoMo1_1 = t[1], FoMo1_2 = t[2]))
  
}


w %>% pivot_longer(-person, names_to = "model", values_to = "weight") %>%
  ggplot(aes(model, weight, fill = factor(person))) +
  geom_col()


