library(tidyverse)

# these functions import data and does some initial processing to get them into 
# our stim and found format. 

# import_data will return a list that contains two dataframes
# d$stimulus contains all of the stimulus information
# d$found tells us which items were selected, and the order in which they were selected. 

# (0, 0) should represent the bottom-left corner of the display.
# as assume the stimulus has width of 1unit, with height depending on the aspect ratio

import_data <- function(dataset, small_test=FALSE) {
  
  # first, find data folder... 
  # we will assume it is in a parent folder of the current working directory
  data_path <- "../"
  
  for (ii in 1:3) {
    
    if ("data" %in% dir(data_path)) {
      data_path <- paste0(data_path, "data/")
      break
    } else {
      data_path <- paste0(data_path, "../")
    }
  }
  
  d <- switch(dataset,
              "clarke2022qjep" = import_clarke2022qjep(data_path, small_test),
              "tagu2022cog"    = import_tagu2022cog(data_path, small_test),
              "tagu2025" = import_tagu2025(data_path, small_test),
              "kristjansson2014plos" = import_kristjansson2014plos(data_path, small_test),
              "hughes2024rsos" = import_hughes2024rsos(data_path, small_test),
              "bhat2025" = import_bhat2025(data_path, small_test),
              "unknown dataset")
  
  return(list(
    dataset = dataset,
    stim = d$stim, 
    found = d$found
  ))
}

get_train_test_split <- function(d) {
  
  test_train_split <- d$found %>% 
    group_by(person, condition) %>% 
    summarise(n = length(unique(trial_p)), .groups = "drop") %>%
    mutate(split =  ceiling((n/2)), .keep = "unused")
  
  training <- list(
    found = d$found %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% 
      filter(trial_p <= split) %>% select(-split),
    stim  = d$stim  %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% 
      filter(trial_p <= split) %>% select(-split))
  
  testing <- list(
    found = d$found %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% 
      filter(trial_p >  split) %>% select(-split),
    stim  = d$stim  %>% full_join(test_train_split, 
                                  by = join_by(person, condition)) %>% 
      filter(trial_p >  split)  %>% select(-split))
 
  testing$found <- fix_person_and_trial(testing$found)
  testing$stim <- fix_person_and_trial(testing$stim)
  training$found <- fix_person_and_trial(training$found)
  training$stim <- fix_person_and_trial(training$stim)
  
  return(list(
    training = training, 
    testing = testing))
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

import_tagu2025 <- function(data_path = "../data/", small_test) {
  
  datafile <- paste0(data_path, "tagu2025/DATA_ALL.csv")
  d <- read_csv(datafile) 
  
  d_tmp <- d %>%
    mutate(condition = Bloc,
           trial = TrialsCompleted,
           person = participant,
           targ_type = img,
           id = TargetID,
           x = x,
           y = y,
           found = nbselect) %>%
    select(condition, trial, person, targ_type, id, x, y, found) %>%
    mutate(condition = as.factor(condition),
           id = parse_number(id),
           item_class = parse_number(targ_type))
  
  if (small_test) {
    
    d_tmp <- filter(d_tmp, person < 10, trial < 10)
  }
  
  d <- fix_person_and_trial(d_tmp)
  
  # scale x to (0, 1) and y to (0, a) where a is the aspect ratio
  
  # first subtract the min
  d %>% mutate(x = x - min(x),
               y = y - min(y)) -> d
  
  xmax <- max(d$x)
  
  d %>% mutate(x = x/xmax,
               y = y/xmax) -> d
  
  # extract stimulus data
  d_stim <- d %>% select(person, condition, trial, id, x, y, item_class, trial_p) %>%
    arrange(person, condition, trial, id) 
  
  # extract behavioral data
  d_found <- d %>% filter(found > 0) %>% 
    arrange(person, condition, trial, found) 
  
  return(list(stim = d_stim,
              found = d_found))
  
  
}

import_tagu2022cog <- function(data_path = "../data/", small_test) {
  
  my_spec <- cols(
    condition = col_double(),
    trial = col_double(),
    observer = col_double(),
    targ_type = col_double(),
    id = col_double(),
    x = col_double(),
    y = col_double(),
    found = col_double())
  
  datafile <- paste0(data_path, "tagu2022/tagu_2020_prox_mouse.csv")
  d <- read_csv(datafile,
                col_types = my_spec) 
  
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
           condition = fct_recode(condition, value = "1", control = "2"),
           item_class = if_else(item_class == 1, 2, 1)) -> d
  
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
  d$y = 1-d$y
  d$y = d$y - min(d$y)
  
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

import_clarke2022qjep <- function(data_path = "../data/", small_test) {
  
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
  
  datafile <- paste0(data_path, "clarke2022qjep/clarke2022qjep_collected.csv")
  
  d <- read_csv(datafile, 
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
  d$y = 1-d$y
  d$y = d$y - min(d$y)
  
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

import_hughes2024rsos <- function(data_path = "../data/", small_test){
  
  # import from csv files. 
  # These were computed by the pre-processing script in data/hughes2024
  
  d_stim <- read_csv(paste0(data_path, "hughes2024rsos/hughes2024rsos_stim.csv"),
                     show_col_types = FALSE) 
  d_found <- read_csv(paste0(data_path, "hughes2024rsos/hughes2024rsos_found.csv"),
                      show_col_types = FALSE)
  
  # fix condition labels
  d_found %>% mutate(
    condition = str_remove(condition, "cond_"),
    condition = str_replace(condition, "conj", "conjunction"),
    condition = str_replace(condition, "_", "_scarce"),
    condition = str_replace(condition, "scarceAB", "equal")) %>%
    separate(condition, c("condition",  "scarcity"), "_") %>%
    mutate(condition = as_factor(condition)) -> d_found
  
  d_stim %>% mutate(
    condition = str_remove(condition, "cond_"),
    condition = str_replace(condition, "conj", "conjunction"),
    condition = str_replace(condition, "_", "_scarce"),
    condition = str_replace(condition, "scarceAB", "equal")) %>%
    separate(condition, c("condition",  "scarcity"), "_") %>%
    mutate(condition = as_factor(condition)) -> d_stim
  
  # fix trial_p numbering
  d_stim %>% 
    mutate(second_half = (trial_p > 5)) %>%
    mutate(blk = as.numeric(factor(scarcity)),
           tp = if_else(second_half, 
                        10 + (blk-1)*5 + trial_p,  
                        (blk-1)*5 + trial_p),
           trial_p = tp) %>%
    select(-second_half, -blk, -tp) %>%
    arrange(person, trial)  -> d_stim
  
  d_found %>% 
    mutate(second_half = (trial_p > 5)) %>%
    mutate(blk = as.numeric(factor(scarcity)),
           tp = if_else(second_half, 
                        10 + (blk-1)*5 + trial_p,  
                        (blk-1)*5 + trial_p),
           trial_p = tp) %>%
    select(-second_half, -blk, -tp) %>%
    arrange(person, trial)  -> d_found
  
  
  return(list(stim = d_stim,
              found = d_found))
}

import_kristjansson2014plos <- function(data_path = "../data/", small_test) {
  
  my_spec <- cols(
    condition = col_double(),
    trial = col_double(),
    observer = col_double(),
    targ_type = col_double(),
    id = col_double(),
    x = col_double(),
    y = col_double(),
    found = col_double())
  
  d <- read_csv(paste0(data_path, "kristjansson2014/human_data_2014_prox.csv"), col_types = my_spec) 
  
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

import_bhat2025 <- function(data_path = "../data/", small_test) {
  
  d_stim <- read_csv(paste0(data_path, "bhat2025/bhat2025_stim.csv"),
                     show_col_types = FALSE) 
  d_found <- read_csv(paste0(data_path, "bhat2025/bhat2025_found.csv"),
                      show_col_types = FALSE)
  d_age <- read_csv(paste0(data_path, "bhat2025/bhat2025_age.csv"),
                      show_col_types = FALSE)
  d_gender <- read_csv(paste0(data_path, "bhat2025/bhat2025_gender.csv"),
                      show_col_types = FALSE)
  d_rt <- read_csv(paste0(data_path, "bhat2025/bhat2025_rt.csv"),
                      show_col_types = FALSE)
  
  d_found %>% 
    mutate(condition = as_factor(condition)) -> d_found
  
  d_stim %>% 
    mutate(condition = as_factor(condition)) -> d_stim

    return(list(stim = d_stim,
                found = d_found,
                age = d_age,
                gender = d_gender,
                rt = d_rt))

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

