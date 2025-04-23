library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/fit_model.R")


iter <- 500

gen_quant("hughes2024rsos", fomo_ver = "1.0", mode = "traintest") 



# fit_model("hughes2024rsos", fomo_ver = "1.3", mode = "traintest",  iter = iter) fit

# fit_model("clarke2022qjep", fomo_ver = "1.0", mode = "traintest",  iter = iter) 
# fit_model("clarke2022qjep", fomo_ver = "1.1", mode = "traintest",  iter = iter) 
# fit_model("clarke2022qjep", fomo_ver = "1.2", mode = "traintest",  iter = iter) 
# fit_model("clarke2022qjep", fomo_ver = "1.3", mode = "traintest",  iter = iter) 
# fit_model("clarke2022qjep", fomo_ver = "1.4", mode = "traintest",  iter = iter) 

gen_quant("kristjansson2014plos", fomo_ver = "1.0", mode = "traintest", iter = iter)
fit_model("kristjansson2014plos", fomo_ver = "1.2", mode = "traintest", iter = iter)


fit_model("tagu2022cog", fomo_ver = "1.0", mode = "traintest", iter = iter)
fit_model("tagu2022cog", fomo_ver = "1.2", mode = "traintest", iter = iter)

#fit_model("tagu2022cog", fomo_ver = "1.3", mode = "traintest", iter = iter)
#fit_model("tagu2022cog", fomo_ver = "1.4", mode = "traintest", iter = iter)

