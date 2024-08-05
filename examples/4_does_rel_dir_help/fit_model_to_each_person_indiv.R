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


mod10 <- cmdstan_model("../../models/simple/FoMo1_0.stan", 
                       cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)


mod11 <- cmdstan_model("../../models/simple/FoMo1_1.stan", 
                     cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)


mod12 <- cmdstan_model("../../models/simple/FoMo1_2.stan", 
                     cpp_options = list(stan_threads = TRUE), force_recompile = TRUE)


iter = 500

#for (pp in 1:24) {
  
  pp = 1

  d_one_person <- filter_one_person(d, pp)
  filename <- paste0("scratch/person", pp)
  
  # prep for model 1.0
  d_list <- prep_data_for_stan(d_one_person$found, d_one_person$stim, c("spatial", "item_class"))
  d_list <- add_priors_to_d_list(d_list, modelver = "1.0")
  d_list$n_trials_to_sim <- 1
  
  # model 1.0
  fit <- mod10$sample(data = d_list,
                      chains = 4, parallel_chains = 4, threads = 4,
                      refresh = 100,
                      init = 1,
                      iter_warmup = iter, iter_sampling = iter,
                      sig_figs = 3)
  
  fit$save_object(paste0(filename, "_10.rds"))
  
  # prep for model 1.1
  d_list <- prep_data_for_stan(d_one_person$found, d_one_person$stim, c("spatial", "item_class"), 
                               remove_last_found = TRUE)
  d_list <- add_priors_to_d_list(d_list, modelver = "1.1")
  d_list$n_trials_to_sim <- 1


  # model 1.1
  fit <- mod11$sample(data = d_list,
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 100,
                    init = 1,
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3)

 fit$save_object(paste0(filename, "_11.rds"))

  # model 1.2
  fit <- mod12$sample(data = d_list,
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 100,
                    init = 1,
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3)

  fit$save_object(paste0(filename, "_12.rds"))
  
#}

rm(fit)

#for (pp in 1:24) {

  filename <- paste0("scratch/person", pp)
  
  m10 <- readRDS(paste0(filename, "_10.rds"))
  m11 <- readRDS(paste0(filename, "_11.rds"))
  m12 <- readRDS(paste0(filename, "_12.rds"))

  loo10 <- m10$loo()
  loo11 <- m11$loo()
  loo12 <- m12$loo()

  saveRDS(loo10, paste0(filename, "_loo10.rds"))
  saveRDS(loo11, paste0(filename, "_loo11.rds"))
  saveRDS(loo12, paste0(filename, "loo12.rds"))

#}
  
d_one_person$found <- d_one_person$found %>%
  filter(found != max(found))
  
d_one_person$found$pareto_k <- loo11$diagnostics$pareto_k 
d_one_person$found %>% filter(pareto_k > 1)

d_one_person$found %>% mutate(pareto_k = if_else(pareto_k > 1, "bad", "good")) %>%
  group_by(found, condition) %>%
  summarise(pareto_mess = sum(pareto_k == "bad")) %>%
  ggplot(aes(found, pareto_mess)) + geom_col(aes(fill = condition)) 

trl = 6

ds <- d_one_person$stim %>%
  filter(trial == trl)

df <- d_one_person$found %>%
  filter(trial == trl) %>%
  mutate(pareto_k = if_else(pareto_k > 1, "bad", "good"))

ggplot(data = ds, aes(x, y)) + 
  geom_point(size = 5, aes(colour = factor(item_class), shape = factor(item_class))) +
  ggrepel::geom_label_repel(data = df, aes(label = found), size = 3) + 
  #scale_colour_manual(values = c(18, 15, 3, 4)) + 
  scale_shape_manual(values = c(15, 19, 3, 4)) +
  geom_path(data = df, aes(colour = pareto_k, group = 1)) 

#w <- tibble()

#for (pp in 1:24) {
  
#  filename <- paste0("scratch/person", pp)
  
#  loo11 <- readRDS(paste0(filename, "loo11.rds"))
#  loo12 <- readRDS(paste0(filename, "loo12.rds"))
  
#  t <- loo::loo_model_weights(list("FoMo 1.1" = loo11, "FoMo 1.2" = loo12))
  
#  w <- bind_rows(w, tibble(
#    person = pp, FoMo1_1 = t[1], FoMo1_2 = t[2]))
  
#}


#w %>% pivot_longer(-person, names_to = "model", values_to = "weight") %>%
#  ggplot(aes(model, weight, fill = factor(person))) +
#  geom_col()

#### testing simulated data

n_trials_per_cond <- 10

n_item_class <- 2
n_item_per_class <- 20
item_class_weights = c(0.5, 0.5, 0, 0)
b_stick = 1
b_memory = 0

abs_dir_tuning = list(kappa = rep(20, 4), theta = c(2, 0.5, 1, 0.5))
rho_delta = 15
rho_psi = 5

d_sim <- sim_foraging_multiple_trials(person = 1, 
                                  condition = "test",
                                  n_trials_per_cond = n_trials_per_cond,
                                  n_item_class =  n_item_class, n_item_per_class = n_item_per_class,
                                  item_class_weights = item_class_weights, item_labels = item_labels,
                                  b_stick = b_stick, 
                                  rho_delta = rho_delta, 
                                  rho_psi = rho_psi, 
                                  abs_dir_tuning = abs_dir_tuning,
                                  b_memory = b_memory,
                                  inital_sel_params = inital_sel_params,
                                  init_sel_lambda = init_sel_lambda)

d_list <- prep_data_for_stan(d_sim$found, d_sim$stim, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.0")
d_list$n_trials_to_sim <- 1

iter = 500

fit <- mod10$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 200, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object(paste0(filename, "_10.rds"))

# prep for model 1.1
d_list <- prep_data_for_stan(d_sim$found, d_sim$stim, c("spatial", "item_class"), remove_last_found = TRUE)
d_list <- add_priors_to_d_list(d_list, modelver = "1.1")
d_list$n_trials_to_sim <- 1


# model 1.1
fit <- mod11$sample(data = d_list,
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 100,
                    init = 1,
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3)

fit$save_object(paste0(filename, "_11.rds"))

# model 1.2
fit <- mod12$sample(data = d_list,
                    chains = 4, parallel_chains = 4, threads = 4,
                    refresh = 100,
                    init = 1,
                    iter_warmup = iter, iter_sampling = iter,
                    sig_figs = 3)

fit$save_object(paste0(filename, "_12.rds"))

m10 <- readRDS(paste0(filename, "_10.rds"))
m11 <- readRDS(paste0(filename, "_11.rds"))
#m12 <- readRDS(paste0(filename, "_12.rds"))

loo10 <- m10$loo()
loo11 <- m11$loo()
#loo12 <- m12$loo()

saveRDS(loo10, paste0(filename, "_loo10.rds"))
saveRDS(loo11, paste0(filename, "_loo11.rds"))

d_sim$found <- d_sim$found %>%
  filter(found != max(found))

d_sim$found$pareto_k <- loo11$diagnostics$pareto_k 
d_sim$found %>% filter(pareto_k > 1)

d_sim$found %>%mutate(pareto_k = if_else(pareto_k > 1, "bad", "good")) %>%
  group_by(found) %>%
  summarise(pareto_mess = sum(pareto_k == "bad")) %>%
  ggplot(aes(found, pareto_mess)) + geom_col() 

