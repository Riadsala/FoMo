dtest <- readRDS("1_fit_models/scratch/d_list/tagu2025/test.rds")
dtest  <- add_priors_to_d_list(dtest, modelver = "1.3")
dtest$grid_offset = 0
fomo_ver_str <- "1_4"
mod <- cmdstan_model(paste0("../models/multi_level/FoMo", fomo_ver_str, ".stan"))

m <- readRDS("1_fit_models/scratch/models/tagu2022cog/train1_3.model")
m_test <- mod$generate_quantities(m, data = dtest, seed = 123)
