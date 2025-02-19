library(cmdstanr)
library(tidyverse)

options(mc.cores = 4)
source("../../functions/import_data.R")


d <- import_data("clarke2022qjep")

d$found$x <- if_else(d$found$x < 0.001, 0.001, d$found$x)
d$found$y <- if_else(d$found$y < 0.001, 0.001, d$found$y)

# take just fist found item on each trial
d$found <- filter(d$found, found == 1)

m <- cmdstan_model("../../models/initial_sel/inital_sel_one_person.stan")

sl <- list(
  N = nrow(d$found),
  L = length(unique(d$found$person)),
  Z = d$found$person,
  n_targets = nrow(d$found),
  K = 2,
  item_x = d$found$x, 
  item_y = d$found$y,
  X = as.numeric(d$found$condition)
)

m <- m$sample(sl, chains = 4)

 m$draws(c("a_x", "b_x", "c_x",
           "a_y", "b_y", "c_y"), format = "df") %>%
   sample_frac(0.05) %>%
   select(-.chain, -.iteration) -> post

 
comp_stuff <- function(a_x, b_x, c_x, a_y, b_y, c_y, .draw) {
  
  x <- seq(0.01, 0.99, 0.01)
  y <- x
  
  llhx1 <- dbeta(x, c_x, c_x)
  llhy1 <- dbeta(y, c_y, c_y)
  
  llhx2 <- dbeta(x, a_x, b_x)
  llhy2 <- dbeta(y, a_y, b_y)
  
  return(tibble(x = x, llhx1, llhy1, llhx2, llhy2, .draw))
 
}


p <- pmap_df(post, comp_stuff)

p %>% pivot_longer(-c(x, .draw), 
                   names_to = "comp", values_to = "z") %>%
ggplot((aes(x, z, colour = comp, group = interaction(comp, .draw)))) + 
  geom_path(alpha = 0.2)
