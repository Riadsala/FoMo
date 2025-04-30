library(tidyverse)
library(cmdstanr)

source("../../functions/import_data.R")
source("../../functions/prep_data.R")
source("../../functions/fit_model.R")

iter <- 100

gen_quant("hughes2024rsos", fomo_ver = "1.0", mode = "traintest", iter_genquant = iter)
gen_quant("hughes2024rsos", fomo_ver = "1.2", mode = "traintest", iter_genquant = iter)
gen_quant("hughes2024rsos", fomo_ver = "1.3", mode = "traintest", iter_genquant = iter)

gen_quant("kristjansson2014plos", fomo_ver = "1.0", mode = "traintest", iter_genquant = iter)
gen_quant("kristjansson2014plos", fomo_ver = "1.2", mode = "traintest", iter_genquant = iter)
gen_quant("kristjansson2014plos", fomo_ver = "1.3", mode = "traintest", iter_genquant = iter)

gen_quant("tagu2022cog", fomo_ver = "1.0", mode = "traintest", iter_genquant = iter)
gen_quant("tagu2022cog", fomo_ver = "1.2", mode = "traintest", iter_genquant = iter)
gen_quant("tagu2022cog", fomo_ver = "1.3", mode = "traintest", iter_genquant = iter)

gen_quant("bhat2025", fomo_ver = "1.0", mode = "traintest", iter_genquant = iter)
gen_quant("bhat2025", fomo_ver = "1.2", mode = "traintest", iter_genquant = iter)
