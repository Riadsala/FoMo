library(tidyverse)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")


dir.create("scratch")
dir.create("scratch/d_list/")

datasets <- c("hughes2024rsos", "kristjansson2014plos", "tagu2022cog")


for (i in 1:length(datasets)){
  
  dataset <- datasets[i]
  print(paste0("currently processing ", dataset))
  
  fomo_preprocess(dataset, delta0 = 20)
  
}


