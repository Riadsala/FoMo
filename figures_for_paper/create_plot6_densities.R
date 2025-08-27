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

scale_colour_brewer_d <- function(..., palette = "Dark2") {
  scale_colour_brewer(..., palette = palette )
}

scale_fill_brewer_d <- function(..., palette = "Dark2") {
  scale_fill_brewer(..., palette = palette)
}

options(
  ggplot2.discrete.colour = scale_colour_brewer_d,
  ggplot2.discrete.fill = scale_fill_brewer_d
)

sf <- "../examples/1_fit_models/scratch"

dataset <- c("hughes2024rsos")   
model_ver <- "1_3"

# read in data
d <- import_data(dataset)

m <- readRDS(paste0(sf, "/models/", dataset, "/fit/", model_ver, ".model"))
post <- extract_post(m, d)


plt_coreparams <- plot_model_fixed(post, nrow = 2)
plt_theta <- plot_model_theta(post, nrow = 1)

plt_coreparams / free(plt_theta) 
ggsave("figs/fig6_post_densities.pdf", width = 5, height = 7.5)

# currently not in paper
plt_indiv <- plot_model_random(post)
