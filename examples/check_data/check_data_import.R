# script to check datasets....

library(tidyverse)

source("../../functions/import_data.R")

################################################################################
# first of all, let's check (x, y) but plotting distributions when found = 1
# I feel it is a sensible assumption that these should cluster in the top left
################################################################################




datasets <- c("clarke2022qjep", "tagu2022cog", "hughes2024rsos", "kristjansson2014plos")



##############################################################
# check distribution of intial item selection
##############################################################

get_first_items <- function(ds) {
  
  d <- import_data(ds)
  
  d$found %>% filter(found == 1) %>%
    mutate(dataset = d$dataset) -> df
  
  return(df)
  
}

d <- map_df(datasets, get_first_items)

d %>% ggplot(aes(x, y, colour = condition)) + 
  geom_point(alpha = 0.5) +
  facet_wrap(~dataset, nrow = 2)  +
  coord_fixed() +
  ggtitle("dist. of first selected items")

##############################################################
# check number of rows going into training and test sets
##############################################################

check_train_test <- function(ds) {
  
  d <- import_data(ds)
  tt <- get_train_test_split(d)
  
  tibble(dataset = ds,
         split = rep(c("total", "train", "test"), each = 2),
         data = rep(c("stim", "found"), 3),
         rows = c(nrow(d$stim), nrow(d$found), 
                  nrow(tt$training$stim), nrow(tt$training$found),
                  nrow(tt$testing$stim), nrow(tt$testing$found))) -> dout
  
  return(dout)
  
}

map_df(datasets, check_train_test) %>%
  pivot_wider(names_from = split, values_from = rows) %>% 
  mutate(sum = train + test) %>%
  knitr::kable()


##############################################################
# check number of trials in dataset
##############################################################

import_data("clarke2022qjep")$found %>%
  filter(found == 1) %>%
  group_by(person, trial_p) %>%
  summarise(n = n())

import_data("hughes2024rsos")$found %>%
  filter(found == 1) %>%
  group_by(person, trial_p) %>%
  summarise(n = n())
