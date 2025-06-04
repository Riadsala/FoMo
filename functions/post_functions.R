library(tidybayes)

# some helpful functions for extracting fixed and random effect posterior of our foraging model 

# cl contains the condition labels

################################################################################################
# extract posterior density estimates
################################################################################################

extract_post <- function(m, d) {
  
  # extract ALL parameters and collect into a list
  # to ensure that we have the same draws for each part,
  # extract everything and then filter

  # First of all, are we dealing with a multi-level model?
  vars <- m$metadata()$stan_variables
  multi_level <- if_else(sum(str_detect(vars, "z_u")) == 1, TRUE, FALSE)

  # get condition labels from the data
  cl <- unique(d$stim$condition)
  
  # extract the main fixed effect components
  post_fixed <- extract_post_fixed(m, cl)
  param_names <- unique(post_fixed$param)
  
  # post_fixed now wants to be in wide format
  post_fixed %>%
    pivot_wider(names_from = "param", 
                values_from = "value") -> post_fixed
  
  # extract prior
  prior <- extract_prior(m)
  
  # collect in a list
  post_list <- list(params = param_names,
                    fixed = post_fixed,
                    prior = prior)
  
  # if we have a multi-level model, extract all the variance parameters
  # and per-person fits
  if (multi_level) {
    
    post_random <- extract_post_random(m, cl)
    post_var <- extract_var(m, cl, param_names)
    
    post_list <- append(post_list, list(random = 
                                          post_random, variances = post_var))
  }
  
  # if theta is in model fit, extract
  if (sum(str_detect(vars, "log_theta")) >= 1) {
    
    post_theta <- extrat_post_theta(m, cl)
    post_list <- append(post_list, list(theta = post_theta))
    
    # if required, also get the individual differences for theta
    if (multi_level) {
      
      post_theta_u <- extrat_post_theta_u(m, cl)
      post_list <- append(post_list, list(utheta = post_theta_u))
      
    }
  }
  
  return(post_list)
  
}

extrat_post_theta_u <- function(m, cl) {
  
  post_absdir <- m$draws("u_log_theta", format = "df") %>%
    as_tibble() %>%
    select(-.chain, -.iteration) %>%
    pivot_longer(starts_with("u_log_theta"), names_to = "comp", values_to = "log_theta") %>%
    mutate(theta = exp(log_theta), .keep = "unused") %>%
    separate(comp, c("condition", "person", "comp"), sep = ",") %>%  
    mutate(condition = factor(condition, labels = cl),
           comp = (parse_number(comp)),
           phi = (comp-1) * pi/2,
           comp = factor(comp)) 
  
  return(post_absdir)
  
}

extrat_post_theta <- function(m, cl) {
  
  # get fixed (mean) effects
  m$draws("log_theta", format = "df") %>%
    as_tibble() %>%
    pivot_longer(starts_with("log_theta"), names_to = "comp", values_to = "log_theta") %>%
    mutate(theta = exp(log_theta), .keep = "unused") %>%
    separate(comp, c("condition", "comp"), sep = ",") %>%  
    mutate(condition = factor(condition, labels = cl),
           comp = (parse_number(comp)),
           phi = (comp-1) * pi/2,
           comp = factor(comp)) %>%
    select(-.chain, -.iteration) -> post_absdir
  
  n_directions = length(unique(post_absdir$phi))
  
  # get variances in direction
  if ("sigma_w" %in% names(m$draws(format = "df"))) {
  
  m$draws("sigma_w", format = "df") %>%
    as_tibble() %>%
    pivot_longer(starts_with("sigma_w"), names_to = "comp", values_to = "sigma") %>%
    mutate(
      comp = parse_number(comp),
      condition = (comp-1) %/% n_directions,
      condition = factor(condition, labels = cl),
           comp = comp %% n_directions,
           comp = factor(if_else(comp == 0, n_directions, comp)),
            ) %>%
    select(-.chain, -.iteration) -> post_theta_sigma
  
  post_absdir <- left_join(post_absdir, post_theta_sigma, 
            by = join_by(.draw, condition, comp))
  
  } 
  
  
  return(post_absdir)
  
}

extract_var <- function(m, cl, param_names) {
  
  post_sig <- m$draws("sigma_u", format = "df")
  
  # work out how many params
  # we -3 due to the .chain, .iteration and .draw variables
  
  # u params are mapped to parameters using + u[4*(kk-1)+1])
  
  n_params_per_condition <- (length(post_sig) - 3) / length(cl)
  
  post_sig %>% as_tibble() %>%
    select(-.chain, -.iteration) %>%
    pivot_longer(-.draw) %>%
    mutate(condition = (parse_number(name)+n_params_per_condition-1) %/% n_params_per_condition,
           condition = factor(condition, labels = cl),
           param = parse_number(name) %% n_params_per_condition,
           param = if_else(param == 0, n_params_per_condition, param),
           param = factor(param, labels = param_names)) %>%
    select(-name) -> post_sig
  
  return(post_sig)
  
}

extract_prior <- function(m) {
  
  # first, get list of variables in the model:  
  vars <- m$metadata()$stan_variables
  # all priors should start with "prior_"
  pvars <- vars[str_detect(vars, "prior_")]
  # now get prior distributions from model object
  prior <- m$draws(pvars, format = "df")
  
  return(prior)
}

extract_post_fixed <- function(m, cl) {
  
  # first, get list of variables in the model:
  vars <- m$metadata()$stan_variables
  # all fixed effects should start with "rho_" or "b_"
  pvars <- vars[str_detect(vars, "^rho_|^b_")]
  # now get prior distributions from model object
  post <- m$draws(pvars, format = "df") %>%
    as_tibble() %>%
    select(-.chain, -.iteration) %>%
    pivot_longer(-.draw, names_to = "param") %>%
    mutate(condition = parse_number(param), 
           condition = factor(condition, labels = cl),
           param = str_remove(param, "\\[[0-9]*\\]")) -> post
  
  return(post)
  
}

extract_post_random <- function(m, cl) {
  
  # first, get list of variables in the model:
  vars <- m$metadata()$stan_variables
  # all fixed effects should start with "u"
  pvars <- vars[str_detect(vars, "^u_")]
  # we don't want to deal with any theta parameters here
  pvars <- pvars[!str_detect(pvars, "u_log_theta")]
  
  post <- m$draws(pvars, format = "df") %>%
    as_tibble() %>%
    select(-.chain, -.iteration) %>%
    pivot_longer(-.draw) %>%
    separate(name, into = c("param", "person"), sep = ",") %>%
    mutate(person = parse_number(person),
           condition = parse_number(param),
           param = str_remove(param, "\\[[0-9]*"),
           condition = factor(condition, labels = cl)) %>%
    pivot_wider(names_from = "param", values_from = "value")
  
  return(post)
  
}

################################################################################################
# extracting post predictions
################################################################################################

extract_pred <- function(dataset, fomo_ver, folder, mode = "split") {
  
  ################################################################################################
  # This function computes the item-by-item accuracy from the model's generated quantities{} block
  # and joins it with the empirical data so that we can measure accuracy
  ################################################################################################
  
  ################################################################################################
  # getting set up, etc
  
  # load in data
  d <- import_data(dataset)
  
  # if we aren't in model==all, then we only need the test data
  if (mode != "all") {
    d <- get_train_test_split(d)
    d <- d$testing
  }
  
  # load model csvfile
  # do we always want the hardcoded -1 here? (ie, chain 1 I think)
  m <- read_cmdstan_csv(paste0(folder, dataset, fomo_ver, "-1.csv"))
  
  # convert to tibble and tidy
  as_tibble(m$generated_quantities) %>%
    mutate(.draw = 1:n()) %>%
    pivot_longer(-.draw, names_to = "param") -> genquant
  
  # extract item-to-item predictions
  itemwise <- extract_item_pred(genquant, d)
  
  # extract whole-trial predictions
  trialwise <- extract_trial_pred(genquant, d$stim)
  trialwiseF <- extract_trial_pred(genquant, d$stim, mode = "F")

  # we no longer need genquant or m
  rm(m, genquant)
    
  return(list(itemwise = itemwise, 
              trialwise = trialwise,
              trialwise_firstfixed = trialwiseF,
              dataset <- dataset,
              model_ver <- fomo_ver))
  
}

extract_item_pred <- function(gq, my_data) {
  
  gq %>%
    rename(P = "value") %>%
    filter(str_detect(param, "P")) %>%
    separate(param, into = c("chain", "param"), sep = "\\.") %>%
    separate(param, into = c("param", "row"), sep = "\\[") %>%
    select(-chain, -param) %>%
    mutate(row = parse_number(row)) -> pred
  
  # check if we should remove the last item to be selected from my_data
  if (nrow(pred) != nrow(my_data$found)) {
    
    max_found = max(my_data$found$found)
    my_data$found %>% filter(found < max_found) -> my_data$found
  }
  
  pred <- full_join(my_data$found %>% mutate(row = 1:n()), 
                    pred, 
                    by = join_by(row)) %>%
    select(-any_of(c("row", "item_class", "x", "y", "rt"))) %>%
    mutate(model_correct = (P == id))
  
  return(pred)
  
}

extract_trial_pred <- function(gq, d_stim, mode = "Q") {
  
  # first, we extract Q from the model
  gq %>%
    rename(id = "value") %>%
    filter(str_detect(param, mode)) %>%
    separate(param, into = c("trial", "found"), sep = ",") %>%
    mutate(trial = parse_number(str_extract(trial, "(?<=\\[)\\d*")),
           found = parse_number(found)) -> sim
  
  # now merge with d_stim to obtain x, y, item_class info etc
  sim %>% full_join(d_stim, by = join_by(trial, id)) %>%
    arrange(.draw, person, trial, found) %>%
    select(.draw, person, condition, trial, trial_p, 
           found, id, item_class, x, y) -> sim

  return(sim)
}

### 
# needs reworked once we decide more on what graphs we present
summarise_acc <- function(post, compute_hpdi = FALSE) {
  
  if ("split" %in% names(post$itemwise)) {
    
    post$itemwise %>% 
      filter(found != 1) %>%
      group_by(split, condition, found, .draw, person, trial) -> post$itemwise
    
  } else {
    
    post$itemwise %>% 
      filter(found != 1) %>%
      group_by(condition, found, .draw, person, trial) -> post$itemwise
    
  }
  
  post$itemwise  %>%
    summarise(trial_acc = mean(model_correct), .groups = "drop_last") %>%
    summarise(person_acc = mean(trial_acc), .groups = "drop_last") %>%
    summarise(accuracy = mean(person_acc), .groups = "drop_last") -> post$itemwise
  
  if (compute_hpdi) {
    
    post$itemwise %>% 
      median_hdci(accuracy, .width = c(0.53, 0.97)) -> post$itemwise
    
  }
  
  # add meta data
  post$itemwise %>% mutate(model = post$model_ver,
                 dataset = post$dataset) -> post$itemwiseacc
  
  if ("split" %in% names(post$itemwise)) {
    
    post$itemwise %>%  mutate(split = factor(split, levels = c("training", "testing"))) -> post$itemwise
    
  }
  
  return(post$itemwise)
  
}
