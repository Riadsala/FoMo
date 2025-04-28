# this script takes a previously fitted model (see 1_fit_models) and creates a load of plots
# giving an overview of how well the model fits the data

library(tidyverse)
library(cmdstanr)
library(patchwork)
library(tidybayes)

source("../functions/import_data.R")
source("../functions/compute_summary_stats.R")
source("../functions/post_functions.R")
source("../functions/plot_model.R")


options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

model_ver <- "1_0"
dataset <- "kristjansson2014plos"

# read in data
d <- import_data(dataset)

sf <- "1_fit_models/scratch"
folder <- paste0(sf, "/post/", dataset, "/")

#############################################################################
# plot accuracy
#############################################################################
pred <- readRDS(paste0(folder, "pred_1_0.rds"))

pred$itemwise %>% 
  group_by(condition, found, .draw, person) %>%
  summarise(person_acc = mean(model_correct), .groups = "drop_last") %>%
  summarise(accuracy = mean(person_acc), .groups = "drop_last") -> acc

plt_acc <- plot_model_accuracy(acc)

rm(acc)

#############################################################################
# compare run statistics
#############################################################################

trl_stats <- read_csv(paste0("1_fit_models/scratch/post/", dataset, "/run_statistics.csv")) %>%
  pivot_longer(starts_with("v"), names_to = "model_version", values_to = "predicted")

# reorder factor levels for plotting
trl_stats %>% 
  mutate(statistic = fct_recode(statistic, 
                                `number of runs` = "num_runs", 
                                `max(run length)` = "max_run_length",
                                pao = "mean_pao",
                                `best r`  = "mean_bestr"),
         statistic = fct_relevel(statistic, 
                                 "number of runs", "max(run length)", "pao")) -> trl_stats

ggplot(trl_stats, aes(observed, predicted, colour = condition)) + 
  geom_point() + 
  geom_abline(linetype = 2) + 
  facet_wrap( ~ statistic, nrow = 2, scales = "free") +
  theme(legend.position = "none") -> plt_runs

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

iisv <- read_csv(paste0("1_fit_models/scratch/post/", dataset, "/iisv_statistics.csv")) %>%
  filter(z %in% c("observed", paste0("v", model_ver))) %>%
  mutate(z = if_else(str_detect(z, model_ver), "predicted", "observed")) %>%
  rename(data = "z")

iisv %>% 
  filter(is.finite(d2)) %>%
  group_by(found, data, condition, person) %>%
  summarise(d2 = mean(d2)) %>%
  summarise(d2 = mean(d2)) %>%
  ggplot(aes(found, d2, colour = data)) +
  geom_path(alpha = 0.5) +
  facet_wrap(~condition) + 
  paletteer::scale_colour_paletteer_d("fishualize::Acanthurus_sohal") +
  theme(legend.position = "none") +
  scale_x_continuous("item selection")-> plt_delta

iisv %>% 
  filter(is.finite(psi)) %>%
  ggplot(aes(psi, fill = data)) + 
  geom_histogram(position = position_identity(),
                 breaks = seq(0, 1, 0.2),
                 alpha = 0.5) + 
  paletteer::scale_fill_paletteer_d("fishualize::Acanthurus_sohal") -> plt_psi

iisv %>% 
  filter(is.finite(theta)) %>%
  ggplot(aes(theta, fill = data)) + 
  geom_histogram(position = position_identity(),
                 breaks = seq(-pi, pi, pi/8),
                 alpha = 0.5) + 
  paletteer::scale_fill_paletteer_d("fishualize::Acanthurus_sohal") -> plt_phi

#############################################################################
#assemble plot
#############################################################################
 
plt_top <- plt_acc / plt_runs + plot_layout(heights = c(2,3), guides = "collect")
plt_bot <- plt_delta / (plt_psi + plt_phi) + plot_layout(guides = "collect")

# ggsave("scratch/old_model_stat_runs.pdf", plt_top, width = 5, height = 6)
# ggsave("scratch/old_model_stat_iisv.pdf", plt_bot, width = 5, height = 3)

#############################################################################
# plot posterior densities
#############################################################################

m <- readRDS(paste0("1_fit_models/scratch/models/", dataset, "/fit/", model_ver, ".model"))
post <- extract_post(m, d)
post_plt <- plot_model_fixed(post)

plot_model_theta(post, nrow = 2)
ggsave("theta_fixed.png", width = 4, height = 6)
