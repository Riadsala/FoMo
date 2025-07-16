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


for (dataset in datasets) {
# create output folder
outfolder <- paste0("scratch/post/", dataset)

iisv <- read_csv(paste0(outfolder, "/iisv_statistics.csv")) %>%
  mutate(model_version = if_else(
    is.na(model_version), z, model_version)) %>%
  select(person, condition, model_version, d2) 

# get levy flight statistic
iisv %>% modelr::data_grid(person, condition, model_version) %>%
  pmap_df(get_levy, iisv) %>%
  pivot_wider(names_from = model_version, values_from = alpha) %>%
  pivot_longer(-c(person, condition, observed), 
               names_to = "model_version", values_to = "alpha") %>%
  mutate(alpha = round(alpha, 3)) -> levy


write_csv(levy, paste0(outfolder, "/levy_flight.csv"))
}
