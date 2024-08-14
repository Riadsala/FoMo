
add_priors_to_d_list <- function(dl, modelver="1.1") {
  
  filename <- paste0("priors_model", modelver, ".csv")
  
  # first, remove any priors that have already been added
  dl <- dl[!str_detect(names(dl), "prior_") ]
  
  # add in priors
  priors <- read_csv(paste0("../../models/multi_level/", filename), show_col_types = FALSE) %>%
    pivot_longer(-param, names_to = "stat", values_to = "value") %>%
    mutate(param = paste("prior", stat, param, sep="_")) %>%
    select(-stat)
  
  dl <- append(dl, priors %>% group_by(param) %>% deframe())
  
  return(dl)
  
}

prep_data_for_stan <- function(df, ds, model_components = "spatial", remove_last_found = FALSE) {
  
  # df and ds should match d$found and d$stim, which are output by import_data()
  # model_components tells us which model_components to include
  

  ###################################################
  # first, do some processing that everything requires
  
  # unsure where to put this... 
  # remove distracters, as current model ignores them
  ds %>% 
     filter(item_class %in% c(1, 2)) -> ds
  
  if (remove_last_found) {
    
    df %>% 
      filter(found != max(found)) -> df
    
  
  }
  # 
  # extract stimulus parameters
  n_people <- length(unique(df$person))
  n_trials <- length(unique(df$trial))
  
  # we will be trying to predict the item IDs, so lets save them as Y
  Y = as.numeric(df$id)
  
  # get condition info
  # this is a little more complicated than it looks as we may have some missing data
  d_trl <- ds %>% group_by(person, trial) %>% 
    summarise(condition = unique(condition), 
              n_items = n(),
              .groups = "drop")
  
  # remove trials in which no targets were found
  d_trl <- filter(d_trl, trial %in% (df %>% group_by(trial) %>% 
                                       summarise(n=n(), .groups = "drop"))$trial)
  
  n_targets = unique(d_trl$n_items)
  
  X <- as.numeric(d_trl$condition)
  
  # add  these to list
  d_list <- list(
    N = nrow(df),
    L = n_people,
    K = length(unique(ds$condition)),
    n_trials = n_trials,
    n_targets = n_targets,
    Y = Y,
    X = array(X),
    Z = df$person,
    found_order = df$found)  
  
  rm(d_trl, X)
  
  ###################################################
  if ("cross_line" %in% model_components) {
    

    df %>% group_by(person, condition, trial) %>%
      mutate(sideh = if_else(x<0.5, 1, 2),
             sidev = if_else(y<0.5, 1, 2)) -> df
    
    ds %>% group_by(person, condition, trial) %>%
      mutate(sideh = if_else(x<0.5, 1, 2),
             sidev = if_else(y<0.5, 1, 2)) -> ds
    
    found_sideh = t(array(as.numeric(ds$sideh), dim = c(max(df$found), n_trials)))
    #found_sideh[which(found_sideh==2)] =  -1
    
    found_sidev = t(array(as.numeric(ds$sidev), dim = c(max(df$found), n_trials)))
    #found_sidev[which(found_sidev==2)] =  -1
    
    matchingh = array()
    matchingv = array()
    trl <- 1
    
    for (ii in 2:length(Y)) {
      if (df$found[ii]==1) {trl = trl + 1}
      

        found_classh <- df$sideh[ii-1]
        found_classv <- df$sidev[ii-1]
        
      
      matchingh <- rbind(matchingh, as.numeric(found_sideh[trl,] == found_classh))
      matchingv <- rbind(matchingv, as.numeric(found_sidev[trl,] == found_classv))
      
    }
    
    matchingh[which(matchingh == 0 )] = -1
    matchingh[which(is.na(matchingh))] = 0
    
    matchingv[which(matchingv == 0 )] = -1
    matchingv[which(is.na(matchingv))] = 0
    
    
    d_list <- append(d_list,
                     list(cross_line_h = matchingh,
                          cross_line_v = matchingv))
  }
  
  
  if ("spatial" %in% model_components) {
  
    # We are assuming (x, y) are already on some sensible scale
    
    # pre-compute direction and distance data
    spatial <- compute_inter_item_directions_and_distances(Y, df, ds)
    
    # pre-compute relative direction data
    rel_direction <- compute_inter_sel_direction(Y, df, ds) 
    
    # rescale x and y to be both in the (0, 1) range
    ds %>% mutate(x = as.vector(scales::rescale(x, to = c(0.01, 0.99))),
                  y = as.vector(scales::rescale(y, to = c(0.01, 0.99)))) -> ds
    
    # get x y coords of all items
    ds %>% ungroup() %>% select(person, trial, id, x) %>%  
      pivot_wider(names_from = "id", values_from = "x") %>%
      select(-person, -trial) -> itemX
    
    ds %>% ungroup() %>% select(person, trial, id, y) %>%  
      pivot_wider(names_from = "id", values_from = "y") %>%
      select(-person, -trial) -> itemY
  
    d_list <- append(d_list, list(item_x = itemX,
                                  item_y = itemY,
                                  delta = spatial$distances,
                                  phi = spatial$directions,
                                  psi = rel_direction))

    rm(spatial, rel_direction)
  }

  if ("item_class" %in% model_components) { 
    
    # Used for categorical model_components, i.e., items are of one class or another
    # We want to extract the number of classes, the number of targets per class, 
    # class ID of each item, and whether the ith selected item matches the i-1 selected item
    
    n_item_class <- length(unique(ds$item_class))
    n_item_per_class <- length(unique(ds$id))/n_item_class
    
    item_class = t(array(as.numeric(ds$item_class), dim = c(n_item_per_class*n_item_class, n_trials)))
    item_class[which(item_class==2)] =  -1
    
    # work out which targets match previous target
    matching <- does_item_match_prev_target(Y, df, item_class, n_item_per_class, n_item_class)
    
    d_list <- append(d_list, list(n_classes = n_item_class,
                                  item_class = item_class,
                                  S = matching))
  }
  
  if ("cts" %in% model_components) {
    
    # Used when we have continuous model_components, i.e, the items have model_components vectors that take on a range of values
    # special exception: circular model_components (colour, orientation) are dealt with separately 
    
    # we need to think asbout variable names
    
    item_cols <- t(array(as.numeric(ds$col), dim = c(n_targets, n_trials)))
    
    d_list <- append(d_list, list(item_colour = item_cols))
    
  }
  
  if ("circ" %in% model_components) {
    
    # Used when we have continuous circular model_components, i.e, values on a colour wheel, orientations
    # For now, we extract the distance matrix for such model_components
    
    col_dist <- compute_inter_item_colour_dist(Y, df, ds)
    
    d_list <- append(d_list, list(col_dist = col_dist))
  }
  
  d_list$trial = df$trial
 
  return(d_list)
  
}

compute_inter_item_directions_and_distances <- function(Y, df, ds) {
  
  directions <- array()
  distances <- array()
  
  for (ii in 1:length(Y)) {
    
    if (df$found[ii]==1) {
      
      #starting a new trial, so get relevant data
      trl_dat = filter(ds, 
                       person == df$person[ii],
                       trial == df$trial[ii])
      
      # as first target in trial, distance and direction will be 0
      # these should be NA, but not supported by Stan
      directions <- rbind(directions, rep(0, nrow(trl_dat)))
      distances <- rbind(distances, rep(0, nrow(trl_dat)))
      
    } else {
      
      x <- trl_dat$x - trl_dat$x[Y[ii-1]]
      y <- trl_dat$y - trl_dat$y[Y[ii-1]]
    
      distances <- rbind(distances, sqrt(x^2 + y^2))
      directions <- rbind(directions, atan2(y, x))
    }
  }
  
  directions <-  directions[-1,]
  distances <-  distances[-1,]

  return(list(directions = directions, 
              distances  = distances))
}

compute_inter_sel_direction <- function(Y, df, ds) {
  
  # function to compute the angular difference between the 
  # vector from target i-2 to i-1 and the vector from 
  # target i-1 to i. 
  
  theta <- array()
  
  for (ii in 1:length(Y)) {
    
    #starting a new trial, so get relevant data
    trl_dat = filter(ds, 
                     person == df$person[ii],
                     trial == df$trial[ii])
    
    if (df$found[ii] %in% c(1, 2)) {
      
      # as first or second target in trial, direction will be zero
      phi <- rep(0, nrow(trl_dat))
      
    } else {
      
      d_targ <- df[ii-1,]
      # 
      # if ((d_targ$item_class == 0)) {
      #   phi <- rep(0, nrow(trl_dat))
      #   
      # } else {
        d_prev_targ <- filter(df, 
                              person == df$person[ii],
                              trial == df$trial[ii],
                              found == df$found[ii] - 2)
        
        phi = atan2((d_targ$y -  d_prev_targ$y), (d_targ$x -  d_prev_targ$x)) * 180/pi
        phi = (atan2((trl_dat$y - d_targ$y), (trl_dat$x - d_targ$x)) * 180/pi) - phi 
        phi = pmin(abs((phi %% 360)), abs((-phi %% 360)))
        phi = phi/180
        #phi[ii] = 1
      # }
    }
    theta <- rbind(theta, phi[trl_dat$id])
  }
  
  theta <-  theta[-1,]
  
  return(theta)
  
}

does_item_match_prev_target <- function(Y, df, item_class, n_item_per_class, n_item_class) {

  matching = array()
  trl <- 1

  for (ii in 2:length(Y)) {
    if (df$found[ii]==1) {trl = trl + 1}

    if (df$item_class[ii-1]== 0) {
      found_class <- rep(0, n_item_per_class*n_item_class)

    } else {
      found_class <- item_class[trl, Y[ii-1]]

    }

    matching <- rbind(matching, as.numeric(item_class[trl,] == found_class))

  }

  matching[which(matching == 0 )] = -1
  matching[which(is.na(matching))] = 0

  return(matching)

}

