library(tidyverse)


# This script reads evaluates model accuracy and computes summaries
# allowing us to compare human and model run statistics and inter-
# item selection vectors. 

source("../../functions/import_data.R")
source("../../functions/compute_summary_stats.R")
source("../../functions/post_functions.R")

options(mc.cores = 4, digits = 2)

draws_for_sim <- 2

############################################################################

datasets <- c("hughes2024rsos", "tagu2022cog", "kristjansson2014plos")  #"clarke2022qjep", "hughes2024rsos", "tagu2022cog",

dataset <- datasets[1]
# create output folder
outfolder <- paste0("scratch/post/", dataset)

iisv <- read_csv(paste0(outfolder, "/iisv_statistics.csv")) %>%
  select(person, condition, z, d2)

# get levy flight statistic
iisv %>% modelr::data_grid(person, condition) %>%
  pmap_df(get_levy, iisv) %>%
  mutate(statistic = "levy") %>%
  rename(observed = "alpha") %>%
  bind_rows(rl) -> rl