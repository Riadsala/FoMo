################################################################################
# Simulate foraging data
#
# The main function is sim_foraging_trial(), which simulates a single trial.
# sim_foraging_multiple_trials() simulates all trials for one person × condition.
# sim_foraging_people() simulates full datasets with random effects.
#
# Helper functions for spatial computations are included at the end.
################################################################################


################################################################################
# Simulate one foraging trial
################################################################################
sim_foraging_trial <- function(
    trl = 1,
    sp = list(
      n_item_class = 4,
      n_item_per_class = c(10, 10, 10, 10),
      item_labels = c("a", "b", "d1", "d2")
    ),
    fp = list(b_a = 0, b_s = 0, rho_delta = 1, rho_psi = 0),
    adp = "off",
    isp = "off",
    items = NULL,
    dev_output = FALSE,
    d0 = 20
) {
  
  ##############################################################################
  # Create stimulus for this trial
  ##############################################################################
  d_stim <- get_stimulus(sp, items, trl)
  
  # Item-class weights from b_a
  item_class_weights <- c(plogis(fp$b_a), 1 - plogis(fp$b_a), 0, 0)
  
  # Initialise dataframe of remaining items
  d_remain <- d_stim %>%
    mutate(
      found = -1, delta = 0, phi = 0, psi = 0,
      prox = 0, rel_dir = 0, abs_dir_tuning = 0,
      b = item_class_weights[d_stim$item_class]
    ) %>%
    mutate(W = b / sum(b))
  
  ##############################################################################
  # Initial item selection
  ##############################################################################
  t <- 1
  
  if (isp == "off") {
    d_found <- sample_n(d_remain, 1, weight = W) %>% mutate(found = 1)
  } else {
    print("Initial selection model not yet implemented.")
  }
  
  d_stim$found[d_found$id[1]] <- 1
  d_remain <- filter(d_remain, id != d_found$id)
  
  ##############################################################################
  # Repeated foraging search
  ##############################################################################
  keep_searching <- TRUE
  
  while (keep_searching) {
    
    t <- t + 1
    
    # Whether next item matches previous class
    match_prev <- if_else(
      d_remain$item_class == d_found$item_class[t - 1],
      1, -1
    )
    
    # Compute distances and angular features
    d_remain <- compute_delta_and_phi(
      d_remain, d_found, t, fp, adp, d0
    )
    
    # Update sampling weights
    d_remain <- d_remain %>%
      mutate(
        W = log(b) + log(plogis(fp$b_s * match_prev)),
        W = W + prox + rel_dir
      )
    
    # Add absolute direction tuning if used
    if (class(adp) == "list") {
      d_remain$W <- d_remain$W + d_remain$abs_dir_tuning
    }
    
    # Normalise
    d_remain$W <- d_remain$W - matrixStats::logSumExp(d_remain$W)
    
    # If probabilities valid, select next
    if (is.finite(sum(exp(d_remain$W)))) {
      
      d_found <- d_found %>% add_row(
        sample_n(d_remain, 1, weight = exp(W))
      )
      
      d_found$found[t] <- t
      d_stim$found[d_found$id[t]] <- t
      
      d_remain <- filter(d_remain, id != d_found$id[t])
      
    } else {
      keep_searching <- FALSE
    }
  }
  
  ##############################################################################
  # Output formatting
  ##############################################################################
  if (!dev_output) {
    d_found <- d_found %>%
      select(trial, found, id, item_class, x, y)
  }
  
  d_stim <- d_stim %>% select(-found)
  d_found <- arrange(d_found, found)
  
  return(list(stim = d_stim, found = d_found))
}



################################################################################
# Simulate multiple trials for a single person × condition
################################################################################
sim_foraging_multiple_trials <- function(
    person = 1,
    condition = "control",
    n_trials_per_cond = 10,
    sp = list(
      n_item_class = 4,
      n_item_per_class = c(10, 10, 10, 10),
      item_labels = c("a", "b", "d1", "d2")
    ),
    fp = list(b_a = 0, b_s = 0, rho_delta = 1, rho_psi = 0),
    adp = "off",
    isp = "off",
    items = NULL,
    dev_output = FALSE,
    d0 = 20,
    rel_proximity = FALSE
) {
  
  # Simulate each trial
  trial_ids <- 1:n_trials_per_cond
  
  sims <- map(
    trial_ids,
    sim_foraging_trial,
    fp = fp,
    sp = sp,
    adp = adp,
    dev_output = dev_output
  )
  
  # Combine "found" tables
  df <- trial_ids %>%
    map_df(~ sims[[.x]]$found) %>%
    mutate(
      condition = as_factor(condition),
      person = person
    ) %>%
    relocate(person, condition)
  
  # Combine stimulus tables
  ds <- trial_ids %>%
    map_df(~ sims[[.x]]$stim) %>%
    mutate(
      condition = as_factor(condition),
      person = person
    ) %>%
    relocate(person, condition)
  
  # Fix trial numbering
  df <- fix_trial_index(df)
  ds <- fix_trial_index(ds)
  
  return(list(stim = ds, found = df))
}



################################################################################
# Simulate an entire experiment: all people × all conditions
################################################################################
sim_foraging_people <- function(
    params,
    rel_proximity = FALSE,
    filename = "sim"
) {
  
  # Extract experiment parameters
  n <- params$e$n_people
  n_conditions <- params$e$n_conditions
  cond_labels <- params$e$condition_labels
  n_trials_per_cond <- params$e$n_trials_per_cond
  
  # Build random effects table
  dpeeps <- tibble(
    person = rep(1:n, n_conditions),
    condition = rep(cond_labels, each = n),
    b_a = rep(params$f$b_a, each = n),
    sd_a = rep(params$v$b_a, each = n),
    b_s = rep(params$f$b_s, each = n),
    sd_s = rep(params$v$b_s, each = n),
    rho_delta = rep(params$f$rho_delta, each = n),
    sd_delta = rep(params$v$rho_delta, each = n),
    rho_psi = rep(params$f$rho_psi, each = n),
    sd_psi = rep(params$v$rho_psi, each = n)
  )
  
  # Apply random effects
  dpeeps <- pmap_df(dpeeps, gen_random_fx, params$s) %>%
    arrange(person, condition)
  
  # Run simulation for every person × condition
  sims <- dpeeps %>%
    pmap(
      sim_person_condition,
      stimulus_params = params$s,
      params$a,
      n_trials = n_trials_per_cond,
      .progress = TRUE
    )
  
  # Combine found
  df <- 1:nrow(dpeeps) %>%
    map_df(~ sims[[.x]]$found) %>%
    mutate(
      trial_p = trial,
      trial = paste(as.numeric(person), as.numeric(condition), trial),
      trial = as.numeric(as_factor(trial))
    )
  
  # Combine stimulus
  ds <- 1:nrow(dpeeps) %>%
    map_df(~ sims[[.x]]$stim) %>%
    mutate(
      trial_p = trial,
      trial = paste(as.numeric(person), as.numeric(condition), trial),
      trial = as.numeric(as_factor(trial))
    )
  
  # Final output structure
  out <- list(
    stim = ds,
    found = df,
    name = filename,
    dp = dpeeps,
    params = params
  )
  
  # Create save folders
  if (!dir.exists("scratch")) dir.create("scratch")
  if (!dir.exists("scratch/data")) dir.create("scratch/data")
  
  saveRDS(out, paste0("scratch/data/", filename, ".RDS"))
  
  return(out)
}



################################################################################
# Run simulation for one person × condition combination
################################################################################
sim_person_condition <- function(
    person, condition,
    b_a, b_s, rho_delta, rho_psi,
    stimulus_params,
    adp,
    n_trials
) {
  
  foraging_params <- list(
    b_a = b_a,
    b_s = b_s,
    rho_delta = rho_delta,
    rho_psi = rho_psi
  )
  
  sim_foraging_multiple_trials(
    person = person,
    condition = condition,
    fp = foraging_params,
    sp = stimulus_params,
    adp = adp,
    n_trials_per_cond = n_trials
  )
}



################################################################################
# Generate one row of random effects for one participant × condition
################################################################################
gen_random_fx <- function(
    person, condition,
    b_a, sd_a,
    b_s, sd_s,
    rho_delta, sd_delta,
    rho_psi, sd_psi,
    stimulus_params
) {
  
  tibble(
    person = person,
    condition = condition,
    b_a = rnorm(1, b_a, sd_a),
    b_s = rnorm(1, b_s, sd_s),
    rho_delta = rnorm(1, rho_delta, sd_delta),
    rho_psi = rnorm(1, rho_psi, sd_psi)
  )
}



################################################################################
# Create or reuse stimulus for a trial
################################################################################
get_stimulus <- function(sp, items, trl) {
  
  if (is.null(items)) {
    d_stim <- gen_stimulus(sp)
  } else {
    d_stim <- items
  }
  
  d_stim$found <- -1
  d_stim$trial <- trl
  
  return(d_stim)
}



################################################################################
# Generate a random spatial stimulus set
################################################################################
gen_stimulus <- function(sp) {
  
  n <- sum(sp$n_item_per_class)
  
  d_stim <- tibble(
    id = 1:n,
    item_class = rep(1:sp$n_item_class, sp$n_item_per_class)
  )
  
  # Random 2D locations
  d_stim$x <- runif(n, 0, 1)
  d_stim$y <- runif(n, 0, 1)
  
  # Normalise to unit square
  d_stim <- d_stim %>%
    mutate(
      x = (x - min(x)) / max(x),
      y = (y - min(y)) / max(y)
    )
  
  return(d_stim)
}



################################################################################
# Compute distance and angular features for remaining items
################################################################################
compute_delta_and_phi <- function(dr, df, t, fp, adp, d0) {
  
  dr <- dr %>%
    mutate(
      delta = d0 * sqrt((df$x[t - 1] - x)^2 + (df$y[t - 1] - y)^2),
      phi = atan2(y - df$y[t - 1], x - df$x[t - 1])
    )
  
  if (t > 2) {
    dr$psi <- atan2(df$y[t - 1] - df$y[t - 2],
                    df$x[t - 1] - df$x[t - 2])
  } else {
    dr$psi <- NA
  }
  
  dr %>%
    mutate(
      psi = wrap_angle_pi(psi - phi)) -> dr
  
  dr <- dr %>%
    mutate(
      prox = -fp$rho_delta * delta,
      rel_dir = if_else(is.finite(psi), -fp$rho_psi * psi, 0)
    )
  
  if (class(adp) == "list") {
    dr <- dr %>%
      mutate(abs_dir_tuning =
               compute_all_von_mises(adp$theta, adp$kappa, phi))
  } else {
    dr <- dr %>%
      mutate(abs_dir_tuning = 0)
  }
  
  return(dr)
}



################################################################################
# Fix trial indexing for combined data tables
################################################################################
fix_trial_index <- function(dat) {
  
  if ("person" %in% names(dat)) {
    
    dat <- dat %>%
      mutate(
        trial_p = trial,
        trial = paste(as.numeric(person), as.numeric(condition), trial)
      )
    
  } else {
    
    dat <- dat %>%
      mutate(
        trial_p = trial,
        trial = paste(as.numeric(condition), trial)
      )
  }
  
  dat <- dat %>%
    mutate(
      trial = as_factor(trial),
      trial = as.numeric(trial)
    )
  
  return(dat)
}



################################################################################
# Check parameter length and replicate if needed
################################################################################
check_and_rep_param <- function(p, r) {
  
  if (length(p) == 1) {
    p <- rep(p, r)
  }
  
  return(p)
}



################################################################################
# Compute absolute-direction tuning using 4 von Mises lobes
################################################################################
compute_all_von_mises <- function(theta, kappa, phi) {
  
  z <-
    compute_von_mises(0,      phi, theta[1], kappa[1]) +
    compute_von_mises(pi/2,   phi, theta[2], kappa[2]) +
    compute_von_mises(pi,     phi, theta[3], kappa[3]) +
    compute_von_mises(3*pi/2, phi, theta[4], kappa[4]) +
    1
  
  return(z)
}



################################################################################
# Compute a von Mises component
################################################################################
compute_von_mises <- function(mu, phi, theta, kappa) {
  
  z <- theta * exp(kappa * cos(phi - mu)) /
    (2 * pi * besselI(kappa, 0))
  
  return(z)
}



################################################################################
# Merge two simple simulation objects
################################################################################
merge_two_simple_d <- function(d1, d2, lab) {
  
  d1_stim <- d1$stim
  n_trial1 <- max(d1_stim$trial)
  
  d1$stim <- d1$stim %>% mutate(condition = lab[1])
  d1$found <- d1$found %>% mutate(condition = lab[1])
  
  d2$stim <- d2$stim %>%
    mutate(condition = lab[2], trial = trial + n_trial1)
  
  d2$found <- d2$found %>%
    mutate(condition = lab[2], trial = trial + n_trial1)
  
  d$stim <- bind_rows(d1$stim, d2$stim)
  d$found <- bind_rows(d1$found, d2$found)
  
  return(d)
}



################################################################################
# Wrap an angle into [-pi, pi] and convert to absolute fraction of pi
################################################################################
wrap_angle_pi <- function(angle) {
  
  a <- angle %% (2 * pi)
  a <- ifelse(a > pi, a - 2 * pi, a)
  
  abs(a) / pi
}
