library(tidyverse)
library(cmdstanr)
library(patchwork)
library(tidybayes)

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

model_ver <- "1_3"
dataset <- "clarke2022qjep"

# read in data
d <- import_data(dataset)

sf <- "1_fit_models/scratch"
folder <- paste0(sf, "/post/", dataset, "/")

#############################################################################
# plot model comparison over models
#############################################################################
plt_hughes <- plot_models_accuracy("hughes2024rsos", scratch_folder = sf) + theme_bw() + ggtitle("Hughes et al (2024, RSOS)")
plt_clarke <- plot_models_accuracy("clarke2022qjep", scratch_folder = sf) + theme_bw() + ggtitle("Clarke et al (2022, QJEP)")


plt_hughes / plt_clarke

# scatter plot of person acc by model
v1 <- "1_0"
v2 <- "1_3"

plot_model_accuracy_comparison(c("hughes2024rsos", "tagu2022cog"), v1, v2, scratch_folder = sf) 
ggsave("acc_comp.png", width = 8, height = 8)
