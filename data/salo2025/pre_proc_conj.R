library(tidyverse)
source("../../functions/import_data.R")

d <- read_csv('E1_conjunction_all_targets_complete.csv')

d_found <- d %>%
  filter(collected == TRUE) %>%
  mutate(trial = trial_number + 1) %>%
  select(person = "participant", condition = "trial_type", trial, target_name = "target_name", 
         id = "target_name", found = "collection_order", x = "world_x",
         y = "world_z") %>%
  separate_wider_delim(id, delim = "Table", names = c("cup", "table")) %>%
  mutate(item_class = str_sub(cup, 1, 4),
         item_class = as.numeric(as.factor(item_class)),
         cupno = parse_number(cup),
         table = as.numeric(table),
         condition = case_when(item_class == 1 ~ "Conjunction1",
                               item_class == 2 ~ "Conjunction2",
                               item_class == 3 ~ "Conjunction2",
                               item_class == 4 ~ "Conjunction1"),
         item_class = case_when(item_class == 1 ~ 1,
                                item_class == 2 ~ 1,
                                item_class == 3 ~ 2,
                                item_class == 4 ~ 2))

# sometimes d_found doesn't seem to be continuous
# can fix this manually here but should check
d_found <- d_found %>%
  group_by(person, trial) %>%
  mutate(found = row_number())

d_stim <- d %>%
  mutate(trial = trial_number + 1) %>%
  select(person = "participant", condition = "trial_type", trial, 
         ida= "target_name", idb = "generic_target_name",  x = "world_x",
         y = "world_z") %>%
  unite("id", ida:idb, na.rm = TRUE) %>%
  mutate(target_name = id) %>%
  separate_wider_delim(id, delim = "Table", names = c("cup", "table")) %>%
  mutate(item_class = str_sub(cup, 1, 4),
         item_class = as.numeric(as.factor(item_class)),
         cupno = parse_number(cup),
         table = as.numeric(table),
         condition = case_when(item_class == 1 ~ "Conjunction1",
                               item_class == 2 ~ "Conjunction2",
                               item_class == 3 ~ "Conjunction2",
                               item_class == 4 ~ "Conjunction1"),
         item_class = case_when(item_class == 1 ~ 1,
                                item_class == 2 ~ 1,
                                item_class == 3 ~ 2,
                                item_class == 4 ~ 2)) %>%
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
  select(person, condition, trial, id, found, item_class, x, y) -> d_found

d_stim %>%
  select(person, condition, trial, id, item_class, x, y) -> d_stim

# scale x to (0, 1) and y to (0, a) where a is the aspect ratio
# swapped round to normal?

# first subtract the min
d_found %>% mutate(x = x - min(x),
                   y = y - min(y)) -> d_found

ymax <- max(d_found$y)

d_found %>% mutate(x = x/ymax,
                   y = y/ymax) -> d_found

d_stim %>% mutate(x = x - min(x),
                  y = y - min(y)) -> d_stim

ymax <- max(d_stim$y)

d_stim %>% mutate(x = x/ymax,
                  y = y/ymax) -> d_stim

d_found$x <- round(d_found$x, 3)
d_found$y <- round(d_found$y, 3)

d_stim$x <- round(d_stim$x, 3)
d_stim$y <- round(d_stim$y, 3)

d_found <- fix_person_and_trial(d_found)
d_stim <- fix_person_and_trial(d_stim)



# saving
write_csv(d_found, "salo2025_found.csv")
write_csv(d_stim, "salo2025_stim.csv")
