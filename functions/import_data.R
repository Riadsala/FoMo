library(tidyverse)

# these functions import data and does some initial processing to get them into 
# our stim and found format. 

# import_data will return a list that contains two dataframes
# d$stimulus contains all of the stimulus information
# d$found tells us which items were selected, and the order in which they were selected. 

# (0, 0) should represent the bottom-left corner of the display.
# as assume the stimulus has width of 1unit, with height depending on the aspect ratio

import_data <- function(dataset, small_test=FALSE) {
  
  d <- switch(dataset,
              "clarke2022qjep" = import_clarke2022qjep(small_test),
              "tagu2022cog"    = import_tagu2022cog(small_test),
              "kristjansson2014plos" = import_kristjansson2014plos(small_test),
              "hughes2024rsos" = import_hughes2024rsos(small_test),
              "unknown dataset")
}



get_train_test_split <- function(d) {
  
  test_train_split <- d$found %>% 
    group_by(person, condition) %>% 
    summarise(n = length(unique(trial_p)), .groups = "drop") %>%
    mutate(split =  ceiling((n/2)), .keep = "unused")
  
  training <- list(
    found = d$found %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% filter(trial_p <= split),
    stim  = d$stim  %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% filter(trial_p <= split))
  
  testing <- list(
    found = d$found %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% filter(trial_p >  split),
    stim  = d$stim  %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% filter(trial_p >  split))
  
  
  
  testing$found <- fix_person_and_trial(testing$found)
  testing$stim <- fix_person_and_trial(testing$stim)
  training$found <- fix_person_and_trial(training$found)
  training$stim <- fix_person_and_trial(training$stim)
  
  
  return(list(training = training, testing = testing))
}

filter_one_person <- function(d, pp) {


  d_found_small <- d$found %>%
    filter(person == pp) %>%
    mutate(person = 1) %>%
    group_by(trial) %>%
    mutate(trial = cur_group_id()) %>%
    ungroup()

  d_stim_small <- d$stim %>%
    filter(person == pp) %>%
    mutate(person = 1)  %>%
    group_by(trial) %>%
    mutate(trial = cur_group_id()) %>%
    ungroup()

  return(list(found = d_found_small, stim = d_stim_small))

}

fix_person_and_trial <- function(d) {
  
  # first arrange data so it has all of person 1, then person 2, etc
  d %>% arrange(person, condition, trial) -> d
  
  # make sure the person index goes from 1 to N with no missing people
  d$person <- as_factor(d$person)
  levels(d$person) <- 1:length(levels(d$person))
  d$person <- as.numeric(d$person)
  
  # make sure trial number is unique over people and conditions
  # save the old trial info in trial_p
  d %>% mutate(
    trial_p = trial,
    trial = paste(as.numeric(person), as.numeric(condition), trial),
    trial = as_factor(trial),
    trial = as.numeric(trial)) -> d
  
  return(d)
  
}

import_tagu2022cog <- function(small_test) {
  
  my_spec <- cols(
    condition = col_double(),
    trial = col_double(),
    observer = col_double(),
    targ_type = col_double(),
    id = col_double(),
    x = col_double(),
    y = col_double(),
    found = col_double())
  
  d <- read_csv("../../data/tagu2022/tagu_2020_prox_mouse.csv", col_types = my_spec) 
  
  # remove trials with NAs 
  na_trls <- filter(d, is.na(targ_type)) %>%
    mutate(key = paste(condition, trial, observer))
  
  d %>% mutate(key = paste(condition, trial, observer)) %>%
    filter(!(key %in% na_trls$key)) %>%
    select(-key) %>% # we don't need this?
    select(person = "observer", condition, trial = "trial",  
           id = "id", found = "found", item_class = "targ_type",
           x = "x", y = "y") %>%
    mutate(item_class = item_class + 1,
           condition = as.factor(condition),
           condition = fct_recode(condition, value = "1", control = "2")) -> d
  
  if (small_test) {
     
     d <- filter(d, person < 10, trial < 10)
   }
  
  d <- fix_person_and_trial(d)
  
  # scale x to (0, 1) and y to (0, a) where a is the aspect ratio
   
  # first subtract the min
  d %>% mutate(x = x - min(x),
               y = y - min(y)) -> d
  
  xmax <- max(d$x)
  
  d %>% mutate(x = x/xmax,
               y = y/xmax) -> d
  
  # flip y coordinates so that (0, 0) is the bottom left
  # d$y = 1-d$y
  # d$y = d$y - min(d$y)
  # 
  # extract stimulus data
  d_stim <- d %>% select(person, condition, trial, id, x, y, item_class, trial_p) %>%
    arrange(person, condition, trial) 
  
  # extract behavioral data
  d_found <- d %>% filter(found > 0) %>% 
    arrange(person, condition, trial, found) 
  
  if (small_test) {
    d_found <- filter(d_found, found == 1)
  }
  
  return(list(stim = d_stim,
              found = d_found))

}

import_clarke2022qjep <- function(small_test) {
  
  my_spec <- cols(
    condition = col_double(),
    trial = col_double(),
    observer = col_double(),
    targ_type = col_double(),
    id = col_double(),
    x = col_double(),
    y = col_double(),
    found = col_double(),
    RT = col_double())
  
  d <- read_csv("../../data/clarke2022qjep/clarke2022qjep_collected.csv", 
                col_types = my_spec)
  
  d  %>%
    select(person = "observer", condition, trial = "trial",  
           id = "id", found = "found", item_class = "targ_type",
           x = "x", y = "y") %>%
    mutate(item_class = item_class + 1,
           condition = as.factor(condition),
           condition = fct_recode(condition, feature = "1", conjunction = "2")) -> d
  
  if (small_test) {

    d <- filter(d, person < 6, trial < 10)
  }
  
  d <- fix_person_and_trial(d)
  
  # flip y coordinates so that (0, 0) is the bottom left
  # d$y = 1-d$y
  # d$y = d$y - min(d$y)
  
  # extract stimulus data
  d_stim <- d %>% select(person, condition, trial, id, x, y, item_class, trial_p) %>%
    arrange(person, condition, trial) 
  
  # extract behavioral data
  d_found <- d %>% filter(found > 0) %>% 
    arrange(person, condition, trial, found) 
  
  if (small_test) {
    d_found <- filter(d_found, found == 1)
  }

  return(list(stim = d_stim,
              found = d_found))
}

import_hughes2024rsos <- function(small_test){
  
  # import from csv files. 
  # These were computed by the pre-processing script in data/hughes2024
  
  d_stim <- read_csv("../../data/hughes2024rsos/hughes2024rsos_stim.csv",
                     show_col_types = FALSE) 
  d_found <- read_csv("../../data/hughes2024rsos/hughes2024rsos_found.csv",
                      show_col_types = FALSE)
  
  # fix condition labels
  d_found %>% mutate(
    condition = str_remove(condition, "cond_"),
    condition = str_replace(condition, "conj", "conjunction"),
    condition = str_replace(condition, "_", "_scarce"),
    condition = str_replace(condition, "scarceAB", "equal"),
    condition = as_factor(condition)) -> d_found
  
  d_stim %>% mutate(
    condition = str_remove(condition, "cond_"),
    condition = str_replace(condition, "conj", "conjunction"),
    condition = str_replace(condition, "_", "_scarce"),
    condition = str_replace(condition, "scarceAB", "equal"),
    condition = as_factor(condition)) -> d_stim
  
  return(list(stim = d_stim,
              found = d_found))
}


import_kristjansson2014plos <- function(small_test) {
  
  my_spec <- cols(
    condition = col_double(),
    trial = col_double(),
    observer = col_double(),
    targ_type = col_double(),
    id = col_double(),
    x = col_double(),
    y = col_double(),
    found = col_double())
  
  d <- read_csv("../../data/kristjansson2014/human_data_2014_prox.csv", col_types = my_spec) 
  
  d %>% 
    select(person = "observer", condition, trial = "trial",  
           id = "id", found = "found", item_class = "targ_type",
           x = "x", y = "y") %>%
    mutate(item_class = item_class + 1,
           condition = as.factor(condition),
           condition = fct_recode(condition, feature = "1", conjunction = "2")) -> d
  
  if (small_test) {
    
    d <- filter(d, person < 10, trial < 10)
  }
  
  d <- fix_person_and_trial(d)
  
  # extract stimulus data
  d_stim <- d %>% select(person, condition, trial, id, x, y, item_class, trial_p) %>%
    arrange(person, condition, trial) 
  
  # extract behavioral data
  d_found <- d %>% filter(found > 0) %>% 
    arrange(person, condition, trial, found) 
  
  if (small_test) {
    d_found <- filter(d_found, found == 1)
  }
  
  return(list(stim = d_stim,
              found = d_found))
  
}

# some functions to sanity check our data... should be run after import_dat, but before anything else
check_d_stim <- function(d) {
  # this function does as much error checking etc as we can
  
  # is there a *person* col? Is it numeric? 
  check_col_name_and_type(d, "person", "numeric")
  # do the same for id, x and y
  check_col_name_and_type(d, "id", "numeric")
  check_col_name_and_type(d, "x", "numeric")
  check_col_name_and_type(d, "y", "numeric")
  # is condition a factor?
  check_col_name_and_type(d, "condition", "factor")
  
  # additional checks
  if ("person" %in% names(d)) {
    
    # check that it we have no missing people etc
    n_peep = length(unique(d$person))
    
    if (!(min(d$person)==1 & max(d$person)==n_peep)) {
      print("are there missing people?")
    }
  }
  
  # check that x and y are in the range [0, 1]
  if (sum(d$x<0) > 0 | sum(d$x>1) > 0 | sum(d$y<0) > 0 | sum(d$y>1) > 0) {
    print("spatial x and/or y out of bounds")
  }
}

check_col_name_and_type <- function(d, n, t) {
  
  # check that the column n exists in dataframe d
  if (n %in% names(d)) {
    
    # make sure it is of the correct type
    if (t=="numeric") {
      if (!is.numeric(d$x)) {
        print("x should be numeric")
      }
    }
    
  } else {
    print(paste("can't find column", n))
  }
}

