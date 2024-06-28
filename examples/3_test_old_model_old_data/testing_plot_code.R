library(tidyverse)
library(cmdstanr)

source("../../functions/post_functions.R")
source("../../functions/plot_model.R")
source("../../functions/import_data.R")

d <- import_data("tagu2022cog")

m <- readRDS("tagu_1_0_tmp.rds")
#

post <- extract_post(m, d, multi_level = TRUE)


plot_model_fixed(post, clist = list(fill = "condition", facet_cond = "condition"))
