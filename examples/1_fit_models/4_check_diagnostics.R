library(tidyverse)
library(cmdstanr)
library(posterior)

# this script (or QMD) checks previously fitted models and assesses how well they fit

dataset <- "kristjansson2014plos"
folder <- paste0("scratch/models/", dataset, "/fit/")
mode <- "train"

models <- unlist(dir(folder))
models <- models[str_detect(models, "model")]

d <- tibble()

for (model in models) {
  
  m <- readRDS(paste0(folder, model))
  
  vars <- m$metadata()$stan_variables
  vars <- vars[!(vars %in% c("P", "Q",  "z_w", "u", "log_lik"))]
  vars <- vars[!str_detect(vars, "prior")]
  # this is quite slow if we have many variables
  rhat <- m$summary(variables = vars)$rhat
  
  d <- bind_rows(d, 
                 tibble(model = model, rhat = rhat))

}


d %>% 
  ggplot(aes(rhat)) +
  geom_histogram(binwidth = 0.01) +
  facet_wrap(~model, scales = "free") + 
  ggtitle("distribution of Rhat statistics")
