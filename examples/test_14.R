library(cmdstanr)
library(tidyverse)
options(mc.cores = 4, digits = 2)
source("../functions/prep_data.R")
source("../functions/import_data.R")
source("../functions/post_functions.R")
source("../functions/compute_summary_stats.R")
source("../functions/plot_data.R")

dtest <- readRDS("1_fit_models/scratch/d_list/tagu2022cog/test.rds")
dtest  <- add_priors_to_d_list(dtest, modelver = "1.3")
dtest$grid_offset = c(0, 0)
fomo_ver_str <- "1_4"
dataset <- "tagu2022cog"
d <- import_data(dataset)

mod <- cmdstan_model(paste0("../models/multi_level/FoMo", fomo_ver_str, ".stan"))
m <- readRDS("1_fit_models/scratch/models/tagu2022cog/train1_3.model")
m_test <- mod$generate_quantities(m, data = dtest, seed = 123)

t <- readRDS(paste0("1_fit_models/scratch/models/tagu2022cog/test1_3.model"))

# get all model predictions
pred <- extract_pred(list(training = m, testing = t), d)

pred$dataset <- dataset
pred$model_ver <- "1_3"

pred14 <- extract_pred(list(training = m, testing = m_test), d)

# compute empirical run statistics
rl <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs = mean(n_runs),
            mean_bestr = mean(best_r),
            mean_pao = mean(pao),
            .groups = "drop") %>%
  mutate(z = "observed")


  # get simulation data for model

  pred <- pred$trialwise %>% filter(.draw == 1)
  
  # compute simulated run statistics
  rlp <- get_run_info_over_trials(pred) %>%
    group_by(person, condition) %>%
    summarise(max_run_length = mean(max_run_length),
              num_runs = mean(n_runs),
              mean_bestr = mean(best_r),
              mean_pao = mean(pao),
              .groups = "drop")
  
# bind everything together
  rl %>% bind_rows(rlp %>% mutate(z = ("v1.3" ))) -> rl
  
  # compute simulated run statistics
  pred <- pred14$trialwise %>% filter(.draw == 1)
  
  rlp <- get_run_info_over_trials(pred) %>%
    group_by(person, condition) %>%
    summarise(max_run_length = mean(max_run_length),
              num_runs = mean(n_runs),
              mean_bestr = mean(best_r),
              mean_pao = mean(pao),
              .groups = "drop")
  
  # bind everything together
  rl %>%
    bind_rows(rlp %>% mutate(z = ("v1.4" ))) -> rl
  
rl  %>% filter(is.finite(person)) %>%
  pivot_longer(c(max_run_length, num_runs, mean_bestr, mean_pao),
               names_to = "statistic") %>%
  pivot_wider(names_from = z, values_from = value) %>%
  pivot_longer(c("v1.3", "v1.4"), names_to = "model_version", values_to = "predicted")-> r
    
  
  ggplot(r, aes(observed, predicted, colour = condition)) + 
    geom_point() + 
    geom_abline(linetype = 2) + 
    ggh4x::facet_grid2(statistic ~ model_version, scales = "free", independent = "all")
  
  # compute correlations
  
  comp_r <- function(condition, statistic, model_version, trl_stats) {
    
    t <- filter(trl_stats, 
                {{condition}} == condition,
                {{statistic}} == statistic, 
                {{model_version}} == model_version)
    
    r <- cor.test(t$observed, t$predicted)$estimate
    
    return(tibble(condition = condition,
                  statistic = statistic,
                  model_version = model_version,
                  r = r))
  }
  
  
  r %>% modelr::data_grid(condition, statistic, model_version) -> to_test
  
  pmap_df(to_test, comp_r, trl_stats = r) %>%
    pivot_wider(names_from = model_version, values_from = r) %>%
    knitr::kable()
  
  r %>% filter(statistic == "mean_bestr", condition == "control", observed > 0.9) 
  
  
  filter(d$found, person == 19)

  
  plot_a_trial(d$stim, d$found, 388
               )
  