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
                             get_sim = TRUE, draw_sample_frac = 0.01)
  
  # add the post predictions into the post list
  post$acc <- pred$acc
  post$sim <- pred$sim
  
  # further summarise accuracy data
  post$acc <- summarise_acc(post, compute_hpdi = FALSE) 
  
  # add metadata
  post$dataset <- d$name
  post$model_ver <- paste0("FoMo", model_ver)
  
  return(post)
  
}

get_rl_and_iisv_statistics <- function(d, sim) {
  
  # compute empirical run statistics
  rle <- get_run_info_over_trials(d$found) %>%
    group_by(person, condition) %>%
    summarise(max_run_length = mean(max_run_length), 
              n_runs         = mean(n_runs),     
              .groups = "drop") %>% 
    mutate(data = "observed")
  
  # compute simulated run statistics
  rlp <- get_run_info_over_trials(sim) %>%
    group_by(person, condition) %>%
    summarise(max_run_length = mean(max_run_length),
              n_runs         = mean(n_runs),    
              .groups = "drop") %>% 
    mutate(data = "simulated")
  
  # bind everything together
  rl <- bind_rows(rle %>% mutate(.draw = 1) ,rlp)
  
  # compute empirical run statistics
  iisve <- get_iisv_over_trials(d$found) %>% 
    mutate(data = "observed")
  
  # compute simulated run statistics
  iisvp <- get_iisv_over_trials(sim) %>%
    mutate(data = "simulated")
  
  # bind everything together
  iisv <- bind_rows(iisve %>% mutate(.draw = 1), iisvp) 
  
  return(list(run_statistics = rl, 
              inter_item_statistics = iisv))
  
}

#################################################################################
# compute interesting stuff
#################################################################################

post <- get_post_and_pred_from_saved_model(d, "1_0", "traintest")
stats <- get_rl_and_iisv_statistics(d, post$sim)

#################################################################################
# some plots
#################################################################################

# accuracy
post$acc %>%
  ggplot(aes(found, accuracy)) +
  stat_lineribbon(alpha = 0.5)

# posterior densities
plot_model_fixed(post, gt = d$params)

# compare iisv statistics
plot_model_human_iisv_comparison(stats$inter_item_statistics)

# compare run statistics
