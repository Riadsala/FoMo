library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")


iter <- 500

fit_model("kristjansson2014plos", fomo_ver = "1.0", mode = "traintest",  iter = iter) 
fit_model("kristjansson2014plos", fomo_ver = "1.1", mode = "traintest",  iter = iter) 

fit_model("tagu2022cog", fomo_ver = "1.0", mode = "traintest",  iter = iter) 
fit_model("tagu2022cog", fomo_ver = "1.1", mode = "traintest",  iter = iter) 

fit_model("clarke2022qjep", fomo_ver = "1.0", mode = "traintest",  iter = iter) 
fit_model("clarke2022qjep", fomo_ver = "1.1", mode = "traintest",  iter = iter) 


fit_model("hughes2024rsos", fomo_ver = "1.0", mode = "traintest",  iter = iter) 
fit_model("hughes2024rsos", fomo_ver = "1.1", mode = "traintest",  iter = iter) 

