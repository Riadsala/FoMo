library(tidyverse)
source("../../functions/import_data.R")

d <- read_csv('E1_feature_all_targets_complete.csv')

d_found <- d %>%
  filter(collected == TRUE) %>%
  mutate(trial = trial_number + 1) %>%
  select(person = "participant", condition = "interactable_type", trial, target_name = "target_name", 
         id = "target_name", found = "collection_order", x = "world_x",
         y = "world_z") %>%
  separate_wider_delim(id, delim = "Table", names = c("cup", "table")) %>%
  mutate(cup = parse_number(cup),
         table = as.numeric(table))

d_stim <- d %>%
  mutate(trial = trial_number + 1) %>%
  select(person = "participant", condition = "interactable_type", trial, 
         ida = "target_name", idb = "generic_target_name", x = "world_x",
         y = "world_z") %>%
  unite("id", ida:idb, na.rm = TRUE) %>%
  mutate(target_name = id) %>%
  separate_wider_delim(id, delim = "Table", names = c("cup", "table")) %>%
  mutate(cup = parse_number(cup),
         table = as.numeric(table)) %>%
  arrange(person, trial, table, cup) %>%
  group_by(person, trial) %>%
  mutate(id = row_number())

# need to get the proper IDs into d_found

for (i in 1:nrow(d_found)) {
  
  person = d_found$person[i]
  cond = d_found$condition[i]
  trial = d_found$trial[i]
  target_name = d_found$target_name[i]
  
  for (j in 1:nrow(d_stim)) {
    
    if ((d_stim$person[j] == person) && (d_stim$condition[j] == cond) && (d_stim$trial[j] == trial) && (d_stim$target_name[j] == target_name)) {
      
      d_found$id[i] <- d_stim$id[j] 
      
      break
      
    }
  }
}

d_found %>%
  select(person, condition, trial, id, found, x, y) -> d_found

d_stim %>%
  select(person, condition, trial, id, x, y) -> d_stim

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

d_found$x <- round(d_found$x, 3)
d_found$y <- round(d_found$y, 3)

d_stim$x <- round(d_stim$x, 3)
d_stim$y <- round(d_stim$y, 3)

d_found <- fix_person_and_trial(d_found)
d_stim <- fix_person_and_trial(d_stim)

# saving
write_csv(d_found, "salo2025_found.csv")
write_csv(d_stim, "salo2025_stim.csv")
