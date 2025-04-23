################################################################################
# Simulate foraging data
#
# the main function is sim_foraging_trial() which does what is says
# 
# sim_foraging_people() allows us to simulate a whole dataset 
# - multiple people, conditions and trials
# - different people have their own random simulation parameters
#
# utility functions are found in subfunctions/spatial_features.R and are shared 
# with prep_data.R (hopefully)
#
################################################################################

sim_foraging_trial <- function(trl = 1, 
                               sp = list(
                                 n_item_class = 4, 
                                 n_item_per_class = c(10, 10, 10, 10), 
                                 item_labels = c("a", "b", "d1", "d2")),
                               fp = list(
                                 b_a = 0, b_s = 0, rho_delta = 1, rho_psi = 0),
                               adp = "off", isp = "off",
                               items = NULL,
                               dev_output = FALSE, 
                               d0 = 20)  
{
  
  ##############################################################################
  # trl is the trial number (useful when generating many trials!)
  # sp contains details for simualting stimulus   
  # fp contains our four main foraging parameters:
  # - b_a, 
  # - b_stick
  # - rho_delta
  # - rho_psi
  # adp is either "off" or contains details for abs dir behaviour
  # isp is either "off" or contains details for initial selection
  # items: pass in pre-generated items (not yet supported)
  # dev_output: include simulation parameters in output
  # d0: delta scaling to get rho_delta onto a sensible scale ~ 1
  ##############################################################################
  
  # first create stimulus
  d_stim <- get_stimulus(sp, items, trl)
  
  # create item class weights based on b_a
  item_class_weights <- c(plogis(fp$b_a), 1 - plogis(fp$b_a), 0, 0)
  
  # initialize dataframe for tracking which items remain:
  d_remain <- d_stim %>%
    mutate(found = -1, delta = 0, phi = 0, psi = 0, 
           prox = 0, rel_dir = 0, abs_dir_tuning = 0, b = item_class_weights[d_stim$item_class])  %>%
    mutate(W = b/sum(b))
  
  ##############################################################################
  # initial item selection
  t <- 1 
  
  if (isp == "off") {
    
    d_found <- sample_n(d_remain, 1, weight = W) %>%
      mutate(found = 1)
    
  } else {
    print("not yet implemented - check code")
    # d_remain %>% mutate(
    #   w1 = dbeta(x, inital_sel_params$a1x, inital_sel_params$b1x) * dbeta(y, inital_sel_params$a1y, inital_sel_params$b1y),
    #   w2 = dbeta(x, inital_sel_params$a2x, inital_sel_params$b2x) * dbeta(y, inital_sel_params$a2y, inital_sel_params$b2y),
    #   W = W * (init_sel_lambda * w1 + (init_sel_lambda - 1) * w2),
    #   W = W / sum(W)) -> d_remain
  }
  
  # make a note that we have found this target
  d_stim$found[d_found$id[1]] <- 1
  
  # remove this point from the stimuli
  d_remain <- filter(d_remain, id != d_found$id)
  
  ##############################################################################
  keep_searching <- TRUE
  
  while(keep_searching) {
    
    t <- t + 1 # we want to find the next item
    
    # compare to item classes to previously selected item class:
    match_prev = if_else(d_remain$item_class == d_found$item_class[t-1], 1, -1)
    
    # calculate distances from each item to the previously select item
    # then compute the prox and direction weights for each item
    d_remain <- compute_delta_and_phi(d_remain, d_found, t,
                                      fp, adp, d0)
    
    d_remain %>% 
      mutate(
        W = log(b) + log(plogis(fp$b_s * match_prev)),
        W = W + prox + rel_dir) -> d_remain
    
    # add in absolute direction tuning if we are using it
    if (class(adp) == "list") d_remain$W = d_remain$W + d_remain$abs_dir_tuning
    
    # normalise
    d_remain$W = d_remain$W - matrixStats::logSumExp(d_remain$W)
    
    # if at least some items have > 0 weight, select one!
    # otherwise, quit search
   if (is.finite(sum(exp(d_remain$W)))) {
      
      # sample the next target
      d_found %>% add_row(
        sample_n(d_remain, 1, weight = exp(W))) -> d_found
      
      d_found$found[t] <- t
      d_stim$found[d_found$id[t]] <- t
      
      # remove target from those that remain
      d_remain <- filter(d_remain, id != d_found$id[t])
      
    } else {
      keep_searching <- FALSE
    }
  }
  
  # unless dev_output == TRUE, remove the sim variables
  if (dev_output == FALSE) {
    
    d_found %>% select(trial, found, id, item_class, x, y) -> d_found
    
  }
  
  # finally, tidy up a little
  # add in sim params to d_stim
  d_stim %>% select(-found) -> d_stim
  d_found %>% arrange(found) -> f_found
  
  return(list(stim = d_stim, found = d_found))
}

################################################################################
# wrapper functions for simulation multiple trials and/or people

sim_foraging_multiple_trials <- function(person = 1,
                                         condition = "control",
                                         n_trials_per_cond = 10,
                                         sp = list(
                                           n_item_class = 4, 
                                           n_item_per_class = c(10, 10, 10, 10), 
                                           item_labels = c("a", "b", "d1", "d2")),
                                         fp = list(
                                           b_a = 0, b_s = 0, rho_delta = 1, rho_psi = 0),
                                         adp = "off", isp = "off",
                                         items = NULL,
                                         dev_output = FALSE, 
                                         d0 = 20,
                                         rel_proximity = FALSE) 
{
  
  ##############################################################################
  # person is the person id (useful when generating many people!)
  # condition is a label   
  # fp contains our four main foraging parameters:
  # - b_a, 
  # - b_stick
  # - rho_delta
  # - rho_psi
  # adp is either "off" or contains details for abs dir behaviour
  # isp is either "off" or contains details for initial selection
  # items: pass in pre-generated items (not yet supported)
  # dev_output: include simulation parameters in output
  # d0: delta scaling to get rho_delta onto a sensible scale ~ 1
  ##############################################################################
  
  # Generate a number of trials with identical sim params. 
  trls <- 1:n_trials_per_cond
  
  d <- map(trls, sim_foraging_trial, 
           fp = fp, sp = sp, adp = adp) 
  
  # rearrange the list structure
  df <- 1:n_trials_per_cond %>% map_df(~d[[.x]]$found) %>%
    mutate(condition = condition,
           condition = as_factor(condition),
           person = person) %>%
    relocate(person, condition)
  
  ds <- 1:n_trials_per_cond %>% map_df(~d[[.x]]$stim) %>%
    mutate(condition = condition,
           condition = as_factor(condition),
           person = person) %>%
    relocate(person, condition)
  
  # correctly label trial and trial_p
  ds <- fix_trial_index(ds)
  df <- fix_trial_index(df)
  
  return(list(stim = ds, found = df))
}

sim_foraging_people <- function(params,
                                rel_proximity = FALSE,
                                filename = "sim") {
  
  # extract experiment parameters
  n <- params$e$n_people
  n_conditions <- params$e$n_conditions
  cond_labels <- params$e$condition_labels
  n_trials_per_cond <- params$e$n_trials_per_cond

  # generate random effects
  dpeeps <- tibble(person      = rep(1:n, n_conditions),
                   condition   = rep(cond_labels, each = n),
                   b_a      = rep(params$f$b_a, each = n),
                   sd_a      = rep(params$v$b_a, each = n),
                   b_s     = rep(params$f$b_s, each = n),
                   sd_s  = rep(params$v$b_s, each = n),
                   rho_delta   = rep(params$f$rho_delta, each = n),
                   sd_delta= rep(params$v$rho_delta, each = n),
                   rho_psi     = rep(params$f$rho_psi, each = n),
                   sd_psi  = rep(params$v$rho_psi, each = n),)
  
  dpeeps <- pmap_df(dpeeps, gen_random_fx, params$s) %>%
    arrange(person, condition)
  
  dpeeps %>% pmap(sim_person_condition, 
                  stimulus_params = params$s,
                  n_trials = n_trials_per_cond, .progress = TRUE) -> d
  
  # rearrange list structure
  df <- 1:nrow(dpeeps) %>% map_df(~d[[.x]]$found) %>%
    mutate(trial_p = trial,
           trial = paste(as.numeric(person), as.numeric(condition), trial),
           trial = as_factor(trial),
           trial = as.numeric(trial))
  
  ds <- 1:nrow(dpeeps) %>% map_df(~d[[.x]]$stim) %>%
    mutate(trial_p = trial,
           trial = paste(as.numeric(person), as.numeric(condition), trial),
           trial = as_factor(trial),
           trial = as.numeric(trial))
  
  d <- list(stim = ds, found = df,
            name = filename,
            dp = dpeeps,
            params = params)
  
  # create save folder if it doesn't yet exist
  
  if(!dir.exists("scratch")) {
    dir.create("scratch")
  }
  
  if(!dir.exists("scratch/data")) {
    dir.create("scratch/data")
  }
  
  saveRDS(d, paste0("scratch/data/", d$name, ".RDS"))
  
  return(d)
  
}

sim_person_condition <- function(person, condition, 
                                 b_a, b_s, rho_delta, rho_psi,
                                 stimulus_params,
                                 n_trials) {
  
  foraging_params <- list(b_a = b_a, 
                          b_s = b_s,
                          rho_delta = rho_delta,
                          rho_psi = rho_psi)
  
  d <- sim_foraging_multiple_trials(person, condition, 
                                    fp = foraging_params,
                                    sp = stimulus_params,
                                    n_trials_per_cond = n_trials)
  
  return(d)
  
}

gen_random_fx <- function(person, condition,
                          b_a, sd_a,
                          b_s, sd_s,
                          rho_delta, sd_delta,
                          rho_psi, sd_psi,
                          stimulus_params) 
{
  

  dout <- tibble(person = person, condition = condition, 
                 b_a = rnorm(1, b_a, sd_a),
                 b_s = rnorm(1, b_s, sd_s),
                 rho_delta = rnorm(1, rho_delta, sd_delta),
                 rho_psi = rnorm(1, rho_psi, sd_psi))
  
  return(dout)
}


################################################################################
# helper functions for simulating a trial

get_stimulus <- function(sp, items, trl) {
  
  # set up dataframe for storing things
  # if we have passed in `items`, use that.
  # otherwise generate a new random stimulus
  if (is.null(items)) {
    
    d_stim <- gen_stimulus(sp)
    
  } else {
    
    d_stim = items
    
  }
  
  d_stim$found = -1
  d_stim$trial = trl
  
  return(d_stim)
  
}

gen_stimulus <- function(sp) {
  
  n <- sum(sp$n_item_per_class)
  
  d_stim <- 
    tibble(
      id = 1:n,
      item_class = rep(1:sp$n_item_class, sp$n_item_per_class))
  
  #class_lab = as_factor(rep(item_labels, n_item_per_class)))
  
  # generate (x, y) locations....
  # we should make this more sophisticated
  d_stim$x <- runif(n, 0, 1)
  d_stim$y <- runif(n, 0, 1)
  
  d_stim %>% mutate(x = x - min(x),
                    y = y - min(y),
                    x = x/max(x),
                    y = y/max(y)) -> d_stim
  
  return(d_stim)
}

compute_delta_and_phi <- function(dr, df, t, fp, adp, d0) {
  
  # for each item, compute - 
  # delta: distance from the previously selected 
  # phi: angle from the previously selected item
  # psi: angle from t-2 selected item to previously selected item 
  
  dr %>% mutate(
    delta = d0 *sqrt((df$x[t-1] - x)^2 + (df$y[t-1] - y)^2),
    phi = 180 * atan2((y - df$y[t-1]), (x - df$x[t-1]))/pi) -> dr
  
  if (t > 2)
  {
    dr$psi <- 180 *atan2((df$y[t-1] - df$y[t-2]), (df$x[t-1] - df$x[t-2])) / pi
  } else { 
    dr$psi <- NA
  }
    
  # psi is updated to be difference between psi and phi
  # shift angles to positive
  dr %>% mutate(psi = psi - phi,
                phi = pmin(abs((phi %% 360)), abs((-phi %% 360))),
                psi = pmin(abs((psi %% 360)), abs((-psi %% 360))),
                psi = psi/180) -> dr
  
  # compute proximity of remaining targets from current target
  dr %>% mutate(
    prox = - fp$rho_delta * delta,
    rel_dir = if_else(is.finite(psi), - fp$rho_psi * psi, 0)) -> dr
  
  if (class(adp) == "list") {
    
    dr %>% mutate(
      abs_dir_tuning = compute_all_von_mises(adp$theta, adp$kappa, phi)) -> dr
    
  } else {
    
    dr %>% mutate(abs_dir_tuning = 0) -> dr
    
  } 
    
  return(dr)
  
}

fix_trial_index <- function(dat) {
  
  if ("person" %in% names(dat)) {
    
    dat %>% mutate(
      trial_p = trial,
      trial = paste(as.numeric(person), as.numeric(condition), trial)) -> dat
    
  } else {
    
    dat %>% mutate(
      trial_p = trial,
      trial = paste(as.numeric(condition), trial)) -> dat
  }
  
  dat %>% mutate(
    trial = as_factor(trial),
    trial = as.numeric(trial)) -> dat
  
  return(dat)
  
}

check_and_rep_param <- function(p, r) {
  
  if (length(p) == 1) {
    p <- rep(p, r)
  }
  
  return(p)
}


compute_all_von_mises <- function(theta, kappa, phi) {
  
  # convert from degrees to radians
  phi <- pi*phi/180
  
  # normalise theta weights
 # theta <- theta/(sum(theta) + 1)
  
  z <- compute_von_mises(0,   phi, theta[1], kappa[1]) + 
       compute_von_mises(pi/2,   phi, theta[2], kappa[2]) +
       compute_von_mises(pi,     phi, theta[3], kappa[3]) +
       compute_von_mises(3*pi/2, phi, theta[4], kappa[4]) +
       1
  
  return(z)
  
}

compute_von_mises <- function(x, phi, theta, kappa) {
  
  z <- theta * exp(kappa * cos(phi-x)) / (2*pi*besselI(kappa,0))
  return(z)
  
}

merge_two_simple_d <- function(d1, d2, lab) {
  
  d1_stim <- d1$stim
  n_trial1 <- max(d1_stim$trial)
  
  d1$stim %>% mutate(condition = lab[1]) -> d1$stim
  d1$found %>% mutate(condition = lab[1]) -> d1$found
  
  d2$stim %>% mutate(condition = lab[2],
                     trial = trial + n_trial1) -> d2$stim
  d2$found %>% mutate(condition = lab[2],
                      trial = trial + n_trial1) -> d2$found
  
  d$stim <- bind_rows(d1$stim, d2$stim)
  d$found <- bind_rows(d1$found, d2$found)
  
  return(d)
  
}