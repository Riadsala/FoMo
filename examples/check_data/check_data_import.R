# script to check datasets....

library(tidyverse)

source("../../functions/import_data.R")

################################################################################
# first of all, let's check (x, y) but plotting distributions when found = 1
# I feel it is a sensible assumption that these should cluster in the top left
################################################################################

get_first_items <- function(ds) {
  
  d <- import_data(ds)
  
  d$found %>% filter(found == 1) %>%
    mutate(dataset = d$dataset) -> df
  
  return(df)
  
}

check_train_test <- function(ds) {
  
  d <- import_data(ds)
  tt <- get_train_test_split(d)
  
  tibble(dataset = ds,
         split = rep(c("train", "test"), each = 2),
         data = rep(c("stim", "found"), 2),
         rows = c(nrow(tt$training$stim), nrow(tt$training$found),
                  nrow(tt$testing$stim), nrow(tt$testing$found))) -> dout
  
  return(dout)
  
}


datasets <- c("clarke2022qjep", "tagu2022cog", "hughes2024rsos", "kristjansson2014plos")

d <- map_df(datasets, get_first_items)

d %>% ggplot(aes(x, y, colour = condition)) + 
  geom_point(alpha = 0.5) +
  facet_wrap(~dataset, nrow = 2)  +
  coord_fixed() +
  ggtitle("dist. of first selected items")
