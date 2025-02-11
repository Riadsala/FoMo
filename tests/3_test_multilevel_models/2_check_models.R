library(tidyverse)
library(cmdstanr)
library(loo)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/fit_model.R")

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

d <- readRDS("scratch/data/test_anna.RDS")

dataset <- "test_anna"
model_ver <- "1_0"
mode <- "all"

get_post_and_pred_from_saved_model <- function(d, model_ver, mode) {
  
  if (mode == "all") {
  
    m <- readRDS(paste0("scratch/models/", d$name, mode, model_ver, ".model"))
    post <- extract_post(m, d)
    
  } else {
    
    m <- readRDS(paste0("scratch/models/", d$name, "train", model_ver, ".model"))
    t <- readRDS(paste0("scratch/models/", d$name, "test", model_ver, ".model"))
    
    # get model posteriors
    post <- extract_post(m, d)
    
    # combine training and testing models
    m <- list(training = m, testing = t)
  }
  
  # get post prediction accuracy
  pred <- summarise_postpred(m, d, 
                             get_sim = TRUE, draw_sample_frac = 0.25)
  
  # add the post predictions into the post list
  post$acc <- pred$acc
  post$acc <- summarise_acc(post, compute_hpdi = FALSE) 
  post$sim <- pred$sim
  
  # further summarise accuracy data
  
  # add metadata
  post$dataset <- d$name
  post$model_ver <- paste0("FoMo", model_ver)
  
  return(post)
  
}


post <- get_post_and_pred_from_saved_model(d, "1_0", "traintest")


# compute empirical run statistics
rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$sim) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) %>%
  pivot_wider(names_from = "x", values_from = "max_run_length") -> rl

write_csv(rl, "scratch/run_statistics_1_0.csv")

iisve <- get_iisv_over_trials(d$found) 

# compute simulated run statistics
iisvp <- get_iisv_over_trials(pred$sim %>%
                                # it would be great to remove this line
                                filter(is.finite(x)))

# bind everything together
bind_rows(iisve %>%  mutate(x = "human"),
          iisvp  %>% mutate(x = "model"))  -> iisv

write_csv(iisv, "scratch/iisv_statistics_1_0.csv")

# 1.1
d <- readRDS("scratch/d_1_1.rds")
m <- readRDS("scratch/sim_train_1_1.model")
t <- readRDS("scratch/sim_test_1_1.model")

# train-test accuracy
pred <- summarise_postpred(list(training = m, testing = t), d, 
                           get_sim = FALSE, draw_sample_frac=0.25)

acc <- compute_acc(pred$acc, compute_hpdi = FALSE) %>% mutate(model = paste("FoMo1.1"))

write_csv(acc, "scratch/post_acc_sim_1.1.csv")

# run length statistics
pred <- summarise_postpred(list(training = m, testing = m), d, 
                           get_sim = TRUE, draw_sample_frac = 0.001)

# compute empirical run statistics
rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$sim) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) %>%
  pivot_wider(names_from = "x", values_from = "max_run_length") -> rl

write_csv(rl, "scratch/run_statistics_1_1.csv")

iisve <- get_iisv_over_trials(d$found) 

# compute simulated run statistics
iisvp <- get_iisv_over_trials(pred$sim %>%
                                # it would be great to remove this line
                                filter(is.finite(x)))

# bind everything together
bind_rows(iisve %>%  mutate(x = "human"),
          iisvp  %>% mutate(x = "model"))  -> iisv

write_csv(iisv, "scratch/iisv_statistics_1_1.csv")

#1.2
d <- readRDS("scratch/d_1_2.rds")
m <- readRDS("scratch/sim_train_1_2.model")
t <- readRDS("scratch/sim_test_1_2.model")

# train test accuracy
pred <- summarise_postpred(list(training = m, testing = t), d, 
                           get_sim = FALSE, draw_sample_frac=0.25)

acc <- compute_acc(pred$acc, compute_hpdi = FALSE) %>% mutate(model = paste("FoMo1.2"))

write_csv(acc, "scratch/post_acc_sim_1.2.csv")

# run length statistics
pred <- summarise_postpred(list(training = m, testing = m), d, 
                           get_sim = TRUE, draw_sample_frac = 0.001)

# compute empirical run statistics
rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$sim) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) %>%
  pivot_wider(names_from = "x", values_from = "max_run_length") -> rl

write_csv(rl, "scratch/run_statistics_1_2.csv")

iisve <- get_iisv_over_trials(d$found) 

# compute simulated run statistics
iisvp <- get_iisv_over_trials(pred$sim %>%
                                # it would be great to remove this line
                                filter(is.finite(x)))

# bind everything together
bind_rows(iisve %>%  mutate(x = "human"),
          iisvp  %>% mutate(x = "model"))  -> iisv

write_csv(iisv, "scratch/iisv_statistics_1_2.csv")


###################################
#### Simulating multiple conditions

item_class_weights = list(c(0.5, 0.5, 0, 0), 
                          c(0.7, 0.3, 0, 0))

b_stick = c(0, 2)

rho_delta = c(20, 15)
sd_rho_delta = 5

rho_psi = c(-1, -1)

abs_dir_tuning = list(kappa = rep(10, 4), theta = rep(1, 4))

d_2cond <- sim_foraging_people(n_people = 12,
                         n_conditions = 2,
                         cond_lab = c("A", "B"),
                         n_trials_per_cond = 8,
                         n_item_class = 2, n_item_per_class = 20,
                         item_class_weights, sd_bA = 0.2,
                         b_stick = b_stick, sd_b_stick = 1,
                         rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                         rho_psi = rho_psi, sd_rho_psi = 0.5,
                         abs_dir_tuning = abs_dir_tuning,
                         inital_sel_params = inital_sel_params) 


saveRDS(d_2cond, "scratch/d_2cond.rds")

# model 1.0

m <- fit_model(d_2cond, fomo_ver = "1.0", mode = "traintest",  iter = 500, n_trials_to_sim = 3) 

# model 1.1

m <- fit_model(d_2cond, fomo_ver = "1.1", mode = "traintest",  iter = 500, n_trials_to_sim = 3) 

# model 1.2

m <- fit_model(d_2cond, fomo_ver = "1.2", mode = "traintest",  iter = 500, n_trials_to_sim = 3) 



#########################################################
## model 1.1

rho_delta = c(1,2)
sd_rho_delta = 0.1

d4 <- sim_foraging_people(n_people = 10,
                          n_conditions = 2,
                          cond_lab = c("A", "B"),
                          n_trials_per_cond = 4,
                          n_item_class = 2, n_item_per_class = 10,
                          item_class_weights, sd_bA = 0.2,
                          b_stick = b_stick, sd_b_stick = 1,
                          rho_delta = rho_delta, sd_rho_delta = sd_rho_delta,
                          rho_psi = rho_psi, sd_rho_psi = 0.5,
                          abs_dir_tuning = abs_dir_tuning,
                          inital_sel_params = inital_sel_params,
                          rel_proximity = TRUE) 

d4$found <- fix_person_and_trial(d4$found)
d4$stim <- fix_person_and_trial(d4$stim)

saveRDS(d4, "scratch/d_2cond_relprox.rds")

d_list <- prep_data_for_stan(d4$found, d4$stim, c("spatial", "item_class"))
d_list <- add_priors_to_d_list(d_list, modelver = "1.1")
d_list$n_trials_to_sim <- 1

iter = 200
mod <- cmdstan_model("../../models/multi_level/FoMo1_1.stan")

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object("scratch/multi_level_2cond_1_1_tmp.rds")

# model 1.2

mod <- cmdstan_model("../../models/multi_level/FoMo1_2.stan")

fit <- mod$sample(data = d_list, 
                  chains = 4, parallel_chains = 4, threads = 4,
                  refresh = 10, 
                  iter_warmup = iter, iter_sampling = iter,
                  sig_figs = 3)

fit$save_object("scratch/multi_level_2cond_1_2_tmp.rds")