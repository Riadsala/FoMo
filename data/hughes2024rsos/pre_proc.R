library(tidyverse)
source("../../functions/import_data.R")

found_spec <- cols(
  person = col_character(),
  block = col_character(),
  condition = col_character(),
  trial = col_double(),
  attempt = col_double(),
  id = col_double(),
  found = col_double(),
  score = col_double(),
  item_class = col_character(),
  x = col_double(),
  y = col_double(),
  rt = col_double())

stim_spec <- cols(
  person = col_character(),
  block = col_character(),
  condition = col_character(),
  trial = col_double(),
  attempt = col_double(),
  id = col_double(),
  item_class = col_character(),
  x = col_double(),
  y = col_double())

# should read in all csvs
p_folders <- dir()

d_stim <- tibble()
d_found <- tibble()
d_age <- tibble()
d_gender <- tibble()

for (pp in 1:length(p_folders)) {
  
  p_file_found <- dir(paste0(p_folders[pp]), "_found.csv", full.names = TRUE)
  p_file_stim <- dir(paste0(p_folders[pp]), "_stim.csv", full.names = TRUE)
  
  p_found <- read_csv(p_file_found, col_types = found_spec)
  p_stim <- read_csv(p_file_stim, col_types = stim_spec)
  
  d_found <- bind_rows(d_found, p_found)
  d_stim <- bind_rows(d_stim, p_stim)
  
  # age and gender
  filename <- str_split_i(p_file_found, "/", 5)
  p_age <- tibble(str_split_i(filename, "_", 2))
  p_gender <- tibble(str_split_i(filename,"_", 3))
  
  d_age <- bind_rows(d_age, p_age)
  d_gender <- bind_rows(d_gender, p_gender)
  
}

# recode person to be a number
d_found %>% mutate(person = parse_number(person)) -> d_found
d_stim %>% mutate(person = parse_number(person)) -> d_stim

# take only the highest number attempt (this needs checking with a dataset with some mistakes)
# I think removing any attempt = 6 should also get rid of any cases with found = 0. This may occasionally throw away a genuine trial?
d_found %>%
  group_by(person, condition, trial) %>%
  filter(attempt == max(attempt)) %>%
  filter(attempt != 6) -> d_found

d_stim %>%
  group_by(person, condition, trial) %>%
  filter(attempt == max(attempt)) %>%
  filter(attempt != 6) -> d_stim

# filter out practice trials (for now?)
d_stim %>%
  filter(condition != "cond_prac") -> d_stim

d_found %>%
  filter(condition != "cond_prac") -> d_found

d_found  %>%
  mutate(trial = trial + 1) %>%
  select(person = "person", condition, trial = "trial",  
         id = "id", found = "found", item_class = "item_class",
         x = "x", y = "y", rt = "rt") %>%
  mutate(item_class = as.numeric(factor(item_class)),
         condition = as.factor(condition),
         condition = as.integer(condition),
         condition = as.factor(condition)) -> d_found

d_stim  %>%
  mutate(trial = trial + 1) %>%
  select(person = "person", condition, trial = "trial",  
         id = "id", item_class = "item_class",
         x = "x", y = "y") %>%
  filter(item_class != "dist_class1") %>%
  filter(item_class != "dist_class2") %>%
  mutate(id_seq = rep(1:20)) %>% # this is hard coded and should be fixed
  mutate(item_class = as.numeric(factor(item_class)),
         condition = as.factor(condition),
         condition = as.integer(condition),
         condition = as.factor(condition)) -> d_stim

# need to get sequential id into d_found

for (i in 1:nrow(d_found)) {
  
  person = d_found$person[i]
  cond = d_found$condition[i]
  trial = d_found$trial[i]
  id = d_found$id[i]
  
  for (j in 1:nrow(d_stim)) {
    
    if ((d_stim$person[j] == person) && (d_stim$condition[j] == cond) && (d_stim$trial[j] == trial) && (d_stim$id[j] == id)) {
      
      d_found$id_seq[i] <- d_stim$id_seq[j] 
      
      break
      
    }
  }
}

d_found %>%
  select(-id) %>%
  rename(id = id_seq) %>%
  select(person, condition, trial, id, found, item_class, x, y, rt) -> d_found

d_stim %>%
  select(-id) %>%
  rename(id = id_seq) %>%
  select(person, condition, trial, id, item_class, x, y) -> d_stim

# scale x to (0, 1) and y to (0, a) where a is the aspect ratio

# first subtract the min
d_found %>% mutate(x = x - min(x),
                   y = y - min(y)) -> d_found

xmax <- max(d_found$x)

d_found %>% mutate(x = x/xmax,
                   y = y/xmax) -> d_found

d_stim %>% mutate(x = x - min(x),
                  y = y - min(y)) -> d_stim

xmax <- max(d_stim$x)

d_stim %>% mutate(x = x/xmax,
                  y = y/xmax) -> d_stim


d_found <- fix_person_and_trial(d_found)
d_stim <- fix_person_and_trial(d_stim)

write_csv(d_found, "hughes2024rsos_found.csv")
write_csv(d_stim, "hughes2024rsos_stim.csv")


