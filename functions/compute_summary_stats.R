library(ForagingOrg)
library(poweRlaw)


# Computes some descriptive summary statistics of foraging data

# get_run_info_over_trials(df) - pass in a dataframe, and get out one row per trial
# - number items found per trial
# - number of runs
# - max run length

# get_inter_sel_infor_over_trials(df) - pass in a dataframe, and ger out one row per item selection
# - inter-item selection distance
# - inter-item selection direction 


get_run_info_over_trials <- function(df) {
  
  # Function to compute run statistics (max run length and number of runs) for 
  # a dataframe (d$found) of trials
  
  # If posterior predictions are supplied, .draw will be included as a column.
  # In this case, we compute the run statistics for each .draw x trial
  
  if (".draw" %in% c(names(df))) {
    
    df %>% unite(trial_p, trial_p, .draw, sep = "_") -> df
    draws_present <- TRUE
    
  } else {
    draws_present <- FALSE
  }

  if ("condition" %in% names(df)) {
    
    df %>% group_by(person, condition, trial_p) -> df
    
  } else {
    
    df %>% group_by(person, trial_p) -> df
    
  }
  
 df %>% 
    summarise(.groups = "drop") %>%
    pmap_df(get_run_info, df = df, 
            .progress = TRUE) -> dout
 
  # now split dout$trial up by draws
  if (draws_present) {
   
    dout %>% separate(trial_p, 
                      into = c("trial_p", ".draw"), 
                      "_", convert = TRUE) -> dout
  }
 
  return(dout)
  
}

get_iisv_over_trials <- function(df) {
  
  # Function to compute iisv statistics for a dataframe (d$found) of trials
  
  # Returns delta (distance), theta (direction) and psi (relative direction)
  
  # If posterior predictions are supplied, .draw will be included as a column.
  # In this case, we compute the run statistics for each .draw x trial

  if (".draw" %in% c(names(df))) {
    
   df %>% unite(trial_p, trial_p, .draw, sep = "_") -> df
    draws_present <- TRUE
    
  } else {
    draws_present <- FALSE
  }
    
  d_trials <- df %>% group_by(person, condition, trial_p) %>% 
    summarise(.groups = "drop") 
  
  dout <- pmap_df(d_trials, get_iisv_stats, df,
                   .progress = TRUE) 
  
  if (draws_present) {
    
    dout %>% separate(trial_p, 
                      into = c("trial_p", ".draw"), 
                      "_",
                      convert = TRUE) -> dout
  }
  
  return(dout)
  
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
    trial_p = trl,
    condition = condition,
    n_found = nrow(trl_dat),
    max_run_length = max(rl$lengths), 
    n_runs = length(rl$lengths),
    best_r = best.r(trl_dat$x, trl_dat$y),
    pao = PAO(trl_dat$x, trl_dat$y)))
}    


get_iisv_stats <- function(person, condition, trial_p, d) {
  
  pp = person
  trl = trial_p
  cnd = condition
  
  trl_dat <- filter(d, person == pp, trial_p == trl,  condition == cnd) %>% 
    arrange(found) %>%
    select(-item_class) %>%
    mutate(x0 = lag(x), y0 = lag(y),
           x00 = lag(x, 2), y00 = lag(y, 2),
           d2 = (x-x0)^2 + (y-y0)^2,
           theta = atan2(y-y0, x-x0),
           theta0 = atan2(y0-y00, x0-x00),
           psi = theta - theta0,
           psi = pmin(abs((psi %% 2*pi)), abs((-psi %% 2*pi))),
           psi = psi / pi) %>%
    select(-x0, -y0, -x00, -y00, -theta0)
  
  return(trl_dat)
}

get_levy <- function(person, condition, iisv) {
  
  d <- sqrt(filter(iisv, 
                   person ==!!person, 
                   condition == !!condition,
                   is.finite(d2))$d2)
  
  m <- conpl(d)
  
  return(tibble(person = person,
                condition = condition,
                alpha = estimate_pars(m)$pars))
  
}


