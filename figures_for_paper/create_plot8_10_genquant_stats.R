library(tidyverse)
library(patchwork)

source("../functions/import_data.R")
# source("../functions/compute_summary_stats.R")
# source("../functions/plot_model.R")
# source("../functions/plot_data.R")
# source("../functions/post_functions.R")
# source("../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

#############################################################################
# generated quantities
#############################################################################
datasets <- c( "tagu2022cog", "hughes2024rsos", "clarke2022qjep", "kristjansson2014plos")

rl <- tibble()

for (ds in datasets) {
  
  rl <- bind_rows(rl,
                  read_csv(paste0("../examples/1_fit_models/scratch/post/", ds, "/run_statistics.csv")) %>%
                    mutate(dataset = ds))
}

rl %>% pivot_longer(-c(person, condition, statistic, observed, dataset), 
                    names_to = "model_ver", values_to = "predicted") %>%
  mutate(dataset = fct_relevel(dataset, c("kristjansson2014plos", "tagu2022cog"))) -> rl

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
            b = summary(lm(observed ~ predicted))$coefficients[2,1], 
            .groups = "drop") -> rl_stats

rl %>% filter(model_ver == "f1_3") %>%
  mutate(statistic = factor(statistic, levels = c("num_runs", "max_run_length", "mean_pao","mean_bestr")),
         statistic = fct_recode(statistic, 
                                `num runs` = "num_runs",
                                `max (run length)` = "max_run_length", 
                                PAO = "mean_pao", 
                                `best-r` = "mean_bestr")) %>%
  ggplot(aes(predicted, observed, colour = condition)) +
  geom_abline(linetype = 2) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(group = 1), colour = "black") +
  ggh4x::facet_grid2(statistic ~ dataset, scales = "free", independent = "all") +
  theme_bw() 

ggsave("figs/fig8_run_stats.pdf", width = 12, height = 8)

rl_stats %>% 
  select(-a, -b) %>%
  mutate(statistic = factor(statistic, levels = c("num_runs", "max_run_length", "mean_pao","mean_bestr"))) %>%
  pivot_wider(names_from = "model_ver", values_from = "r") %>%
  knitr::kable()               

# 
# rl  %>%
#   group_by(dataset, statistic, model_ver) %>%
#   summarise(r = cor.test(predicted, observed)$estimate, 
#             a = summary(lm(observed ~ predicted))$coefficients[1,1],
#             b = summary(lm(observed ~ predicted))$coefficients[2,1],.groups = "drop") -> rl_stats



rl %>% mutate(abs_err = abs(observed - predicted), .keep = "unused") %>%
  mutate(first_selection = if_else(str_detect(model_ver, "f"), "fixed", "free"),
         model_ver = str_remove(model_ver, "v|f"),
         model_ver = str_replace(model_ver, "_", "."),
         model_ver = as.numeric(model_ver)) %>%
  group_by(dataset, model_ver, first_selection, statistic) %>%
  summarise(abs_err = median(abs_err)) %>%
  mutate(first_selection = factor(first_selection),
         first_selection =  fct_relevel(first_selection, "free"),
         statistic = factor(statistic),
         statistic = fct_relevel(statistic, "num_runs", "max_run_length")) %>%
  # pivot_wider(names_from = "first_selection", values_from = "tot_err") %>%
  filter(model_ver == "1.3") %>%
  # mutate(percent = 100*(fixed-free)/(free)) %>%
  ggplot(aes(dataset, abs_err, fill = first_selection)) + 
  geom_col(position = position_dodge(), alpha = 0.75) +
  theme_bw() +
  facet_wrap(~statistic, scales = "free_y") +
  scale_x_discrete(guide = guide_axis(angle = 45))

ggsave("figs/fig10_fixed_v_free.pdf", width = 6, height = 5)
  
