library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)

source("../../functions/plot_model.R")

# this script (or QMD) checks previously fitted models and assesses how well they fit

# which dataset and model version do you want?
dataset <- "kristjansson2014plos"
model_ver <- "1_0"

# some set up
folder <- paste0("scratch/models/", dataset, "/fit/")
mode <- "train"
model <- paste0(model_ver, ".model")
m <- readRDS(paste0(folder, model))

vars <- m$metadata()$stan_variables
vars <- vars[!(vars %in% c("P", "Q",  "z_w", "u", "log_lik", "lp__"))]
vars <- vars[!str_detect(vars, "prior")]
# this is quite slow if we have many variables
dsummary <- m$summary(variables = vars)


############################################
# DIAGNOSTIC PLOTS #########################
############################################

# check NA rhats are unused values in L:

print(filter(dsummary, is.na(rhat))$variable)

# check distribution of Rhat statistics
dsummary %>% filter(is.finite(rhat)) %>%
  ggplot(aes(rhat)) +
  geom_histogram() +
  ggtitle("distribution of Rhat statistics")

# traceplots
plot_traceplots(m)

# pairplot
mcmc_pairs(m$draws(), pars = c("b_a[1]", "b_s[1]", "rho_delta[1]", "rho_psi[1]"))
