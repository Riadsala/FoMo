library(tidyverse)
library(patchwork)
library(circular)

source("../functions/import_data.R")
source("../functions/prep_data.R")
source("../functions/compute_summary_stats.R")
source("../functions/plot_model.R")
source("../functions/plot_data.R")
source("../functions/post_functions.R")
source("../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

sf <- "../examples/1_fit_models/scratch"

dataset <- c("hughes2024rsos")   
v2 <- "1_3"

