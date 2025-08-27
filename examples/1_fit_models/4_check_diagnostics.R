library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

# this script (or QMD) checks previously fitted models and assesses how well they fit

dataset <- "kristjansson2014plos"
folder <- paste0("scratch/models/", dataset, "/fit/")
mode <- "train"

models <- unlist(dir(folder))
models <- models[str_detect(models, "model")]



model <- models[1]

m <- readRDS(paste0(folder, model))

vars <- m$metadata()$stan_variables
vars <- vars[!(vars %in% c("P", "Q",  "z_w", "u", "log_lik", "lp__"))]
vars <- vars[!str_detect(vars, "prior")]
# this is quite slow if we have many variables
dsummary <- m$summary(variables = vars)

# check NA rhats are unused values in L:
print(filter(dsummary, is.na(rhat))$variable)

dsummary %>% filter(is.finite(rhat)) %>%
  ggplot(aes(rhat)) +
  geom_histogram() +
  ggtitle("distribution of Rhat statistics")

#######################
# trankplot?

mcmc_trace(m$draws(), regex_pars = "b_a")
mcmc_trace(m$draws(), regex_pars = "b_s")


mcmc_pairs(m$draws(), pars = c("b_a[1]", "b_s[1]", "rho_delta[1]", "rho_psi[1]"))
