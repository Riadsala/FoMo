# Simulate foraging data

# the main function is sim_foraging_trial() which does what is says.
# sim_foraging_people() allows us to simulate a whole dataset (multiple people, conditions and trials)

sim_foraging_people <- function(n_people = 4,
                                n_conditions = 2,
                                cond_lab = c("A", "B"),
                                n_trials_per_cond = 3,
                                n_item_class = 3, n_item_per_class = 2,
                                item_class_weights, sd_bA,
                                b_stick = 0, sd_b_stick = 1,
                                rho_delta = 10, sd_rho_delta = 2,
                                rho_psi = 0, sd_rho_psi = 5,
                                abs_dir_tuning = abs_dir_tuning,
                                inital_sel_params,
                                rel_proximity = FALSE,
                                filename = "sim") {
  
  ## if some params have been specified as a constant, replicate over conditions.
  # b_stick <- check_and_rep_param(b_stick, n_conditions)
  # rho_delta <- check_and_rep_param(rho_delta, n_conditions)
  # sig_theta <- check_and_rep_param(sig_theta, n_conditions)
  # n_trials_per_cond <- check_and_rep_param(n_trials_per_cond, n_conditions)
  
  # generate random effects
  dpeeps <- tibble(person = rep(1:n_people, n_conditions),
                   condition  = rep(cond_lab, each = n_people),
                   mu_cw = rep(item_class_weights, each = n_people),
                   sd_bA = sd_bA,
                   b_stick = rep(b_stick, each = n_people),
                   sd_b_stick = sd_b_stick,
                   rho_delta = rep(rho_delta, each = n_people),
                   sd_rho_delta = sd_rho_delta,
                   rho_psi = rep(rho_psi, each = n_people),
                   sd_rho_psi = sd_rho_psi)
  
  dpeeps <- pmap_df(dpeeps, gen_random_fx)
  
  d <- pmap(dpeeps, sim_foraging_multiple_trials,
            n_trials_per_cond = n_trials_per_cond,
            n_item_class = n_item_class, n_item_per_class = n_item_per_class,
            abs_dir_tuning = abs_dir_tuning,
            inital_sel_params = inital_sel_params,
            item_labels = item_labels, b_memory = 0,
            rel_proximity = rel_proximity,
            .progress = TRUE)
  
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
            dp = dpeeps)
  
  # create save folder if it doesn't yet exist
  if(!dir.exists("scratch/data")) {
    dir.create("scratch/data")
  }
  
  saveRDS(d, paste0("scratch/data/", d$name, ".RDS"))
  
  return(d)
  
}

gen_random_fx <- function(person, condition,
                          mu_cw, sd_bA,
                          b_stick, sd_b_stick,
                          rho_delta, sd_rho_delta,
                          rho_psi, sd_rho_psi,
                          inital_sel_params) 
{
  
  mu_cw[1] <- mu_cw[1] + rnorm(1, 0, sd_bA)
  mu_cw[1] <- if_else(mu_cw[1]<0, 0, mu_cw[1])
  
  init_sel_lambda <- runif(1)
  
  dout <- tibble(person, condition,
                 item_class_weights = list(mu_cw),
                 b_stick = rnorm(1, b_stick, sd_b_stick),
                 rho_delta = rnorm(1, rho_delta, sd_rho_delta),
                 rho_psi = rnorm(1, rho_psi, sd_rho_psi),
                 init_sel_lambda = init_sel_lambda)
  
  return(dout)
}

sim_foraging_multiple_trials <- function(person = 1,
                                         condition = "control", n_trials_per_cond = 10,
                                         n_item_class =  n_item_class, n_item_per_class = n_item_per_class,
                                         item_class_weights = item_class_weights, item_labels =  item_labels,
                                         b_stick = b_stick, 
                                         rho_delta = rho_delta, 
                                         rho_psi = rho_psi, 
                                         abs_dir_tuning = abs_tuning,
                                         b_memory,
                                         inital_sel_params,
                                         init_sel_lambda = 0.25,
                                         rel_proximity = FALSE) 
{
  
  # Generate a number of trials with identical sim params. 
  
  trls <- 1:n_trials_per_cond
  
  d <- map(trls, sim_foraging_trial, 
           n_item_class =  n_item_class, 
           n_item_per_class = n_item_per_class,
           item_class_weights = item_class_weights, 
           item_labels =  item_labels,
           b_stick = b_stick, 
           rho_delta = rho_delta, 
           rho_psi = rho_psi, 
           abs_dir_tuning = abs_dir_tuning,
           b_memory = b_memory,
           inital_sel_params = inital_sel_params,
           init_sel_lambda = init_sel_lambda,
           rel_proximity = rel_proximity) 
  
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


sim_foraging_trial <- function(trl = 1, 
                               n_item_class = 4, 
                               n_item_per_class = c(10, 10, 10, 10), 
                               item_class_weights = c(0.5, 0.5, 0, 0),
                               item_labels = c("A", "B", "d1", "d2"),
                               b_stick = 0, # stick weights
                               rho_delta = 0, 
                               rho_psi = 0,
                               abs_dir_tuning = list(kappa = rep(0, 4), theta = rep(1, 4)),
                               omi_dir = 0,
                               b_memory = 0,
                               inital_sel_params,
                               init_sel_lambda,
                               items = NULL,
                               dev_output = FALSE,
                               rel_proximity = FALSE)  
{
  
  # n_class is the number of different target classes
  # n_per_type is the number of target's per class
  n_item_per_class <- check_and_rep_param(n_item_per_class, n_item_class)
  
  # item_class_weights is the salience score for target type
  item_class_weights <- check_and_rep_param(item_class_weights, n_item_class)
  # normalise
  item_class_weights <- item_class_weights / sum(item_class_weights)
  
  # b_stick is the stick v switch preference
  # rho_delta and sig_theta define the spatial bias
  # trl is the trial number
  
  # calculate total number of items
  n <- sum(n_item_per_class) 
  
  d_stim <- get_stimulus(items, n, n_item_class, n_item_per_class)
  d_stim$trial <- trl
  
  # pick a first point at random, 
  # b is based only on the item_class_weights
  
  d_remain <- d_stim %>%
    mutate(found = -1, delta = 0,  phi = 0, psi = 0, 
           prox = 0, rel_dir = 0, abs_dir_tuning = 0, b = item_class_weights[d_stim$item_class])  %>%
    mutate(W = b/sum(b),
           Wprev = 0,
           Wnew = 0)
  
  # pick a first point at random
  t <- 1 # t is for the "t-th target selection
  
  # initial item selection
  # d_remain %>% mutate(
  #   w1 = dbeta(x, inital_sel_params$a1x, inital_sel_params$b1x) * dbeta(y, inital_sel_params$a1y, inital_sel_params$b1y),
  #   w2 = dbeta(x, inital_sel_params$a2x, inital_sel_params$b2x) * dbeta(y, inital_sel_params$a2y, inital_sel_params$b2y),
  #   W = W * (init_sel_lambda * w1 + (init_sel_lambda - 1) * w2),
  #   W = W / sum(W)) -> d_remain
   
  d_found <- sample_n(d_remain, 1, weight = W) %>%
    mutate(found = 1)
  
  # make a note that we have found this target
  d_stim$found[d_found$id[1]] <- 1
  
  # remove this point from the stimuli
  d_remain <- filter(d_remain, id != d_found$id)
  
  # Now decide if we are stopping already:
  # this could be replaced with a more sophisticated stopping rule
  keep_searching <- TRUE
  
  while(keep_searching) {
    
    t <- t + 1 # we want to find the next item
    
    # compare to previous target
    prev_targ <- d_found$item_class[t-1]
    match_prev = if_else(d_remain$item_class == prev_targ, 1, -1)
    
    # calculate distances from each item to the previously select item
    # then compute the prox and direction weights for each item
    d_remain <- compute_delta_and_phi(d_remain, d_found, t,
                                      rho_delta, abs_dir_tuning, rho_psi,
                                      rel_proximity)
    
    d_remain %>% 
      mutate(
        Wprev = W,
        W = b * plogis(b_stick * match_prev),
        W = W * prox * rel_dir * abs_dir_tuning) -> d_remain
    
    # new weights should be a weighted sum of old weights and current weights
    d_remain %>% 
      mutate(Wnew = W + b_memory * Wprev) -> d_remain
    
    # if at least some items have > 0 weight, select one!
    # otherwise, quit search
    
    if (sum(d_remain$Wnew) > 0) {
      
      # normalise selection weights
      d_remain$Wnew = d_remain$Wnew / sum(d_remain$Wnew) 
      
      # sample the next target
      d_found %>% add_row(
        sample_n(d_remain, 1, weight = Wnew)) -> d_found
      
      d_found$found[t] <- t
      d_stim$found[d_found$id[t]] <- t
      
      d_remain <- filter(d_remain, id != d_found$id[t])
      
      # remove target from those that remain
      
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

check_and_rep_param <- function(p, r) {
  
  if (length(p) == 1) {
    p <- rep(p, r)
  }
  
  return(p)
}

compute_delta_and_phi <- function(dr, df, t, pt, adt, rdt, rel_proximity) {
  
  
  # for each item, compute - 
  # delta: distance from the previously selected 
  # phi: angle from the previously selected item
  # psi: angle from t-2 selected item to previously selected item 
  
  dr %>% mutate(
    delta = sqrt((df$x[t-1] - x)^2 + (df$y[t-1] - y)^2),
    phi = 180 * atan2((y - df$y[t-1]), (x - df$x[t-1]))/pi) -> dr
  
  if (t > 2) {
    
    dr$psi = 180 *atan2((df$y[t-1] - df$y[t-2]), (df$x[t-1] - df$x[t-2])) / pi
    
  } else {
    
    dr$psi = NA
  }
  
  # psi is updated to be difference between psi and phi
  # shift angles to positive
  dr %>% mutate(psi = psi - phi,
                phi = pmin(abs((phi %% 360)), abs((-phi %% 360))),
                psi = pmin(abs((psi %% 360)), abs((-psi %% 360))),
                psi = psi/180) -> dr
  
  # convert to rel_proximity if toggle is on
  if (rel_proximity & nrow(dr)>0 ) {
    
    min_delta <- min(dr$delta)
    dr$delta <- dr$delta / min_delta
  }
  
  # compute proximity of remaining targets from current target
  dr %>% mutate(
    prox = exp(-pt * delta),
    rel_dir = if_else(is.finite(psi), exp(-rdt * psi), 1),
    abs_dir_tuning = compute_all_von_mises(adt$theta, adt$kappa, phi)) -> dr
  
  return(dr)
  
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

get_stimulus <- function(items, n, n_item_class, n_item_per_class) {
  
  # set up dataframe for storing things
  # if we have passed in `items`, use that.
  # otherwise generate a new random stimulus
  if (is.null(items)) {
    
    d_stim <- gen_stimulus(n, n_item_class, n_item_per_class)

  } else {
    
    d_stim = items
    
  }
  
  d_stim$found = -1
  
  return(d_stim)
  
}

gen_stimulus <- function(n, n_item_class, n_item_per_class) {
  
  d_stim <- 
    tibble(
      id = 1:n,
      item_class = rep(1:n_item_class, n_item_per_class))
  
  #class_lab = as_factor(rep(item_labels, n_item_per_class)))
  
  # generate (x, y) locations....
  # we should make this more sophisticated
  d_stim$x <- runif(n, 0, 1)
  d_stim$y <- runif(n, 0, 1)
  
  d_stim %>% mutate(  x = x - min(x),
                      y = y - min(y),
                      x = x/max(x),
                      y = y/max(y)) -> d_stim
  
  return(d_stim)
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