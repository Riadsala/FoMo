library(tidyverse)
library(cmdstanr)
library(loo)

source("../../functions/sim_foraging_data.R")
source("../../functions/post_functions.R")
source("../../functions/prep_data.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/fit_model.R")

#############################################################################
# some set up
#############################################################################

dataset <- "ac_test"
d <- readRDS(paste0("scratch/data/", dataset, ".RDS"))
model_ver <- "1_3"
folder <- paste0("scratch/post/", dataset, "/")

#############################################################################
# plot posterior densities
#############################################################################
m <- readRDS(paste0("scratch/models/", dataset, "/train", model_ver, ".model"))
post <- extract_post(m,d)

post_plt <- plot_model_fixed(post, gt = d$params)


#############################################################################
# plot model accuracy
#############################################################################

# plot model accuracy
plot_models_accuracy(dataset)

# scatter plot of person acc by model
v1 <- "1_0"
v2 <- "1_5"

plot_model_accuracy_comparison(dataset, v1, v2)

