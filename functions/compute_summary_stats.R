# Computes some descriptive summary statistics of foraging data


# get_run_info_over_trials(df) - pass in a dataframe, and get out one row per trial
# - number items found per trial
# - number of runs
# - max run length

# get_inter_sel_infor_over_trials(df) - pass in a dataframe, and ger out one row per item selection
# - inter-item selection distance
# - inter-item selection direction 


get_run_info_over_trials <- function(df) {
  
  df %>%
    group_by(person, condition, trial_p) %>% 
    summarise(.groups = "drop") %>%
    pmap_df(get_run_info, df = df, 
            .progress = TRUE) -> dout
  
  return(dout)
  
}

get_inter_sel_info_over_trials <- function(df) {
  
  d_trials <- d$found %>% group_by(person, condition, trial_p) %>% 
    summarise(.groups = "drop") 
  
  d_out <- pmap_df(d_trials, get_inter_targ_stats, df,
                   .progress = TRUE) 
  
  return(d_out)
  
}

get_run_info <- function(person, condition, trial_p, df) {
  
  # calculate some simple run statistics for a trial
  pp = person
  cnd = condition
  trl = trial_p
  
  trl_dat <- filter(df, person == pp, trial_p == trl,  condition == cnd) %>% 
    arrange(found)
  
  rl <- rle(trl_dat$item_class)
  
  return(list(
    person = person,
    trial = trl,
    condition = condition,
    n_found = nrow(trl_dat),
    max_run_length = max(rl$lengths), 
    n_runs = length(rl$lengths)))
}    

get_inter_targ_stats <- function(person, condition, trial_p, d) {
  
  pp = person
  trl = trial_p
  cnd = condition
  
  trl_dat <- filter(d, person == pp, trial_p == trl,  condition == cnd) %>% 
    arrange(found) %>%
    select(-item_class) %>%
    mutate(x0 = lag(x), 
           y0 = lag(y),
           d = sqrt((x-x0)^2 + (y-y0)^2),
           theta = atan2(y-y0, x-x0)) %>%
    filter(found > 1) %>%
    select(-x0, -y0)
  
  return(trl_dat)
}



