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
  if (sum(str_detect(vars, "log_theta")) > 1) {
    
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
  
  post_absdir <- m$draws("log_theta", format = "df") %>%
    as_tibble() %>%
    pivot_longer(starts_with("log_theta"), names_to = "comp", values_to = "log_theta") %>%
    mutate(theta = exp(log_theta), .keep = "unused") %>%
    separate(comp, c("condition", "comp"), sep = ",") %>%  
    mutate(condition = factor(condition, labels = cl),
           comp = (parse_number(comp)),
           phi = (comp-1) * pi/2,
           comp = factor(comp))
  
  return(post_absdir)
  
}

extract_var <- function(m, cl, param_names) {
  
  post_sig <- m$draws("sigma_u", format = "df")
  
  # work out how many params
  # we -3 due to the .chain, .iteration and .draw variables
  n_params_per_condition <- (length(post_sig) - 3) / length(cl)
  
  post_sig %>% as_tibble() %>%
    select(-.chain, -.iteration) %>%
    pivot_longer(-.draw) %>%
    mutate(condition = parse_number(name) %% length(cl),
           condition = if_else(condition == 0, length(cl), condition),
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

extract_pred <- function(m, d) {
  
  # This function computes the item-by-item accuracy from the model's generated quantities{} block
  # and joins it with the empirical data so that we can measure accuracy
  
  # first, deal with m being either one model or a list
  if (unique(class(m)=="list")) {
    
    # m has been input as a list, presumably:
    # i) a model to to training data
    # ii) test set predictions
    mode <- "train_test"
    mtr <- m$training
    mte <- m$testing
    
  } else {
    
    mode <- "all"
    mtr <- m
    training <- d
    
  }
  
  # get list of variables in the model:
  vars <- mtr$metadata()$stan_variables
  
  # get item-selection predictions
  # sometimes we might want to extract W (loglik) 
  # not done at present as naming inconsistent
  pvars <- vars[str_detect(vars, "^P")]
  
  # do the trial level simulations exist?
  get_sim = ("Q" %in% vars)
  
  # determine if m is a multi-level model or not
  multi_level <- ("z_u" %in% vars)
  
  if (mode == "train_test") {
    
    # get the training/test data split
    dtt <- get_train_test_split(d)
    training <- dtt$training
    testing <- dtt$testing
    rm(dtt)
    
    # get training set predictions
    pred_tr <- extract_item_pred(mtr, training, pvars) %>% mutate(split = "training")
    # now we also need to get the test set predictions
    pred_te <- extract_item_pred(mte, testing, pvars) %>% mutate(split = "testing")
    
    pred <- bind_rows(pred_tr,pred_te)
    
    rm(pred_tr, pred_te)
    
  } else {
    
    # get training set predictions
    pred <- extract_item_pred(mtr, d, pvars)
    
  }
  
  if (get_sim) {
    
    sim <- extract_trial_pred(mtr, training$stim)
    
    # define output list
    list_out <- list(itemwise = pred, trialwise = sim)
    
  } else {
    
    # define output list
    list_out <- list(acc = pred)
  }
  
  return(list_out)
  
}

extract_item_pred <- function(my_model, my_data, pv, chains = 1) {
  
  my_model$draws(pv, format = "df") %>% 
    filter(.chain == chains) %>%
    # sample_frac(sample_frac) %>%
    as_tibble() %>%
    select(-.chain, -.iteration) %>%
    pivot_longer(-.draw) %>%
    mutate(n = parse_number(name),
           name = str_remove(name, "\\[\\d*\\]")) %>%
    pivot_wider(names_from = "name", values_from = "value") -> pred
  
  pred <- full_join(my_data$found %>% mutate(n = 1:n()), 
                    pred, by = join_by(n)) %>%
    select(-any_of(c("n", "item_class", "x", "y", "rt"))) %>%
    mutate(model_correct = (P == id))
  
}

extract_trial_pred <- function(mod, d_stim, mode, chains = 1) {
  
  # do not call directly..
  # this is run by summarise_postpred()
  
  # first, we extract Q from the model
  sim <- mod$draws("Q", format = "df")  %>%
    filter(.chain == chains) %>%
    as_tibble() %>%
    select(-.chain, -.iteration) %>%
    pivot_longer(-.draw, values_to = "id") %>%
    separate(name, 
             c("trial", "found"), 
             sep = ",", convert = TRUE) %>%
    mutate(trial = parse_number(trial), 
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
