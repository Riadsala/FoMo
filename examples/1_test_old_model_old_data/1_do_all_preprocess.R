library(tidyverse)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")


dir.create("scratch")
dir.create("scratch/d_list/")


datasets <- c("hughes2024rsos", "clarke2022qjep")

dataset <- "hughes2024rsos"



fomo_preprocess(dataset)

