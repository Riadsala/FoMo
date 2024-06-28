library(tidyverse)
library(cmdstanr)

source("../../functions/post_functions.R")
source("../../functions/import_data.R")

d <- import_data("tagu2022cog")

m <- readRDS("tagu_1_0_tmp.rds")
#

post <- summarise_postpred(m, d)
