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

model_ver <- "1_3"
dataset <- "hughes2024rsos"

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

plot_model_accuracy_comparison("hughes2024rsos", v1, v2, scratch_folder = sf) -> plt_acc_comp

#############################################################################
# improvement throughout a trial?
#############################################################################

acc <- tibble()
  
  folder <- paste0(sf, "/post/", dataset, "/")
  
  acc1 <- readRDS(paste0(folder, "pred_", v1, ".rds"))$itemwise %>%
    mutate(version = v1) %>%
    filter(found > 1, found < 40) %>%
    group_by(version, condition, found, .draw) %>%
    summarise(accuracy = mean(model_correct),
              .groups = "drop_last") %>%
    median_hdci(accuracy)
  
  acc2 <- readRDS(paste0(folder, "pred_", v2, ".rds"))$itemwise %>%
    mutate(version = v2) %>%
    filter(found > 1, found < 40) %>%
    group_by(version, condition, found, .draw) %>%
    summarise(accuracy = mean(model_correct),
              .groups = "drop_last") %>%
    median_hdci(accuracy)
  
  acc <-bind_rows(acc1, acc2)

rm(acc1, acc2)

acc %>% 
  mutate(version = str_replace(version, "_", "."),
         version = factor(version, levels = c("1.3", "1.0"))) %>%
  ggplot(aes(found, ymin = .lower, ymax = .upper,  fill = version)) +
                 geom_lineribbon(alpha = 0.5) +
  facet_wrap(~condition) +
  scale_x_continuous("item", limits = c(1, 20), expand = c(0,0)) +
  scale_y_continuous("accuracy") -> plt_found

plt_acc_comp / plt_found
ggsave("acc_comp.pdf", width = 6, height = 5)


#############################################################################
# generated quantities
#############################################################################

rl <- tibble()

for (ds in datasets) {
  
  rl <- bind_rows(rl,
                  read_csv(paste0("1_fit_models/scratch/post/", ds, "/run_statistics.csv")) %>%
                    mutate(dataset = ds))
}

rl %>% pivot_longer(-c(person, condition, statistic, observed, dataset), 
                    names_to = "model_ver", values_to = "predicted") -> rl

comp_corr <- function(df) {
  
  r <- cor.test(df$value, df$observed)$estimate
  
  return(tibble(
    dataset = df$dataset, 
    statistic = df$statistic,
    r = r
  ))
}

rl  %>%
  group_by(dataset, statistic, model_ver, condition) %>%
  summarise(r = cor.test(predicted, observed)$estimate, 
            a = summary(lm(observed ~ predicted))$coefficients[1,1],
            b = summary(lm(observed ~ predicted))$coefficients[2,1],.groups = "drop") -> rl_stats

ds <- "tagu2022cog"

rl %>% filter(model_ver == "v1_3") %>%
  mutate(statistic = factor(statistic, levels = c("num_runs", "max_run_length", "mean_pao","mean_bestr")),
         statistic = fct_recode(statistic, `num runs` = "num_runs", `max (run length)` = "max_run_length", PAO = "mean_pao", `best-r` = "mean_bestr")) %>%
  ggplot(aes(predicted, observed, colour = condition)) +
  geom_abline(linetype = 2) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggh4x::facet_grid2(statistic ~ dataset, scales = "free", independent = "all") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("run_stats.pdf", width = 8, height = 9)

rl_stats %>% 
  select(-a, -b) %>%
  mutate(statistic = factor(statistic, levels = c("num_runs", "max_run_length", "pao","mean_bestr"))) %>%
  pivot_wider(names_from = "model_ver", values_from = "r") %>%
  knitr::kable()               
