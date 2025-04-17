library(tidyverse)
library(cmdstanr)
library(patchwork)
library(tidybayes)

source("../functions/import_data.R")
source("../functions/prep_data.R")
source("../functions/compute_summary_stats.R")
source("../functions/plot_model.R")
source("../functions/plot_data.R")
source("../functions/post_functions.R")
source("../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

model_ver <- "1_3"
dataset <- "clarke2022qjep"

# read in data
d <- import_data(dataset)

sf <- "1_fit_models/scratch"
folder <- paste0(sf, "/post/", dataset, "/")

#############################################################################
# plot model comparison over models
#############################################################################
plt_hughes <- plot_models_accuracy("hughes2024rsos", scratch_folder = sf) + theme_bw() + ggtitle("Hughes et al (2024, RSOS)")
plt_clarke <- plot_models_accuracy("clarke2022qjep", scratch_folder = sf) + theme_bw() + ggtitle("Clarke et al (2022, QJEP)")


plt_hughes / plt_clarke

# scatter plot of person acc by model
v1 <- "1_0"
v2 <- "1_3"

plot_model_accuracy_comparison(c("hughes2024rsos", "clarke2022qjep"), v1, v2, scratch_folder = sf) 
ggsave("acc_comp.png", width = 8, height = 8)

#############################################################################
# plot accuracy
#############################################################################
acc <- read_csv(paste0(folder, "acc_train1_0.csv"))
plot_model_accuracy(acc)
ggsave("clarke1_0_acc.png", width = 6, height = 4)
rm(acc)

#############################################################################
# plot posterior densities
#############################################################################
dataset <- "hughes2024rsos"

m <- readRDS(paste0("1_fit_models/scratch/models/", dataset, "/train", model_ver, ".model"))
post <- extract_post(m, d)
post_plt <- plot_model_fixed(post)

plot_model_theta(post, nrow = 2)
ggsave("theta_fixed.png", width = 4, height = 6)

plot_model_theta(post, per_person = TRUE, nrow = 10)
ggsave("test.png", width = 10, height = 20)

#############################################################################
# create plot
#############################################################################
acc_plt / post_plt

#############################################################################
# compare run statistics
#############################################################################

trl_stats <- read_csv(paste0("1_fit_models/scratch/post/", dataset, "/run_statistics.csv")) %>%
  pivot_longer(starts_with("v"), names_to = "model_version", values_to = "predicted")

ggplot(r, aes(observed, predicted, colour = condition)) + 
  geom_point() + 
  geom_abline(linetype = 2) + 
  ggh4x::facet_grid2(statistic ~ model_version, scales = "free", independent = "all")

# compute correlations

comp_r <- function(condition, statistic, model_version, trl_stats) {
  
  t <- filter(trl_stats, 
              {{condition}} == condition,
              {{statistic}} == statistic, 
              {{model_version}} == model_version)
  
  r <- cor.test(t$observed, t$predicted)$estimate
  
  return(tibble(condition = condition,
                statistic = statistic,
                model_version = model_version,
                r = r))
}


trl_stats %>% modelr::data_grid(condition, statistic, model_version) -> to_test

pmap_df(to_test, comp_r, trl_stats = trl_stats) %>%
  pivot_wider(names_from = model_version, values_from = r) %>%
  knitr::kable()

#############################################################################
# compute & compare iisv statistics
#############################################################################


