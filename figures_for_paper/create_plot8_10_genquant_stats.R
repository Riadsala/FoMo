library(tidyverse)
library(patchwork)

source("../functions/import_data.R")
# source("../functions/compute_summary_stats.R")
# source("../functions/plot_model.R")
# source("../functions/plot_data.R")
# source("../functions/post_functions.R")
# source("../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

scale_colour_brewer_d <- function(..., palette = "Dark2") {
  scale_colour_brewer(..., palette = palette )
}

scale_fill_brewer_d <- function(..., palette = "Dark2") {
  scale_fill_brewer(..., palette = palette)
}

options(
  ggplot2.discrete.colour = scale_colour_brewer_d,
  ggplot2.discrete.fill = scale_fill_brewer_d
)

#############################################################################
# generated quantities
#############################################################################

datasets <- c("hughes2024rsos", "kristjansson2014plos", "tagu2022cog", "clarke2022qjep")

rl <- tibble()

for (ds in datasets) {
  
  rl <- bind_rows(rl,
                  read_csv(paste0("../examples/1_fit_models/scratch/post/", ds, "/run_statistics.csv")) %>%
                    mutate(dataset = ds))
  
  # add in levy
  read_csv(paste0("../examples/1_fit_models/scratch/post/", ds, "/levy_flight.csv")) %>%
    rename(predicted = "alpha") %>%
    mutate(dataset = ds, statistic = "levy") %>%
    pivot_wider(names_from = model_version, values_from = "predicted") %>%
    bind_rows(rl) -> rl
}



rl %>%
  filter(is.na(.draw)| .draw == "1") %>%
  select(-.draw) %>% 
  pivot_longer(-c(person, condition, statistic, observed, dataset), 
                     names_to = "model_ver", values_to = "predicted") %>%
  group_by(person, condition, dataset, model_ver, statistic) %>%
  summarise(
    observed = observed, 
    predicted = mean(predicted)) -> rl

rl %>% filter(!str_detect(model_ver, "k")) -> rl

comp_corr <- function(df) {
  
  r <- cor.test(df$value, df$observed)$estimate
  
  return(tibble(
    dataset = df$dataset, 
    statistic = df$statistic,
    r = r
  ))
}

rl <- rl %>%
  filter(model_ver != ".draw")

rl  %>%
  group_by(dataset, statistic, model_ver, condition) %>%
  summarise(r = cor.test(predicted, observed)$estimate, 
            a = summary(lm(observed ~ predicted))$coefficients[1,1],
            b = summary(lm(observed ~ predicted))$coefficients[2,1], 
            .groups = "drop") -> rl_stats


rl %>% filter(model_ver == "v1_3") %>% ungroup() %>%
  mutate(statistic = factor(statistic, levels = c("num_runs", "max_run_length", "mean_pao","mean_bestr", "levy", "mean_int")),
         statistic = fct_recode(statistic, 
                                `num runs` = "num_runs",
                                `max (run length)` = "max_run_length", 
                                PAO = "mean_pao", 
                                `best-r` = "mean_bestr",
                                `mean intersections` = "mean_int"),
         condition = fct_relevel(condition, "conjunction", "feature", "control", "value")) %>%
  ggplot(aes(predicted, observed, colour = condition)) +
  geom_abline(linetype = 2) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F, aes(group = 1), colour = "black") +
  ggh4x::facet_grid2(dataset ~statistic, scales = "free", independent = "all") +
  theme_bw() 

ggsave("figs/fig8_run_stats.pdf", width = 12, height = 8)

rl_stats %>% 
  select(-a, -b) %>%
  mutate(statistic = factor(statistic, levels = c("num_runs", "max_run_length", "mean_pao","mean_bestr", "mean_int"))) %>%
  pivot_wider(names_from = "model_ver", values_from = "r") %>%
  knitr::kable()               

#### FIG 10

source("../functions/plot_model.R")
  
# having a go at scaling things

rl_z <- rl %>%
  mutate(abs_err = abs(observed - predicted), .keep = "unused") %>%
  ungroup() %>%
  #reframe(scale_abs_err = scale(abs_err), .by = c(dataset, statistic)) %>%
  mutate(model_ver = rl$model_ver) %>%
  mutate(first_selection = if_else(str_detect(model_ver, "f"), "fixed", "free"),
         model_ver = str_remove(model_ver, "v|f"),
         model_ver = str_replace(model_ver, "_", "."),
         model_ver = as.numeric(model_ver)) %>%
  group_by(dataset, model_ver, first_selection) %>%
  summarise(abs_err = mean(abs_err)) 

rl_z %>% 
  filter(model_ver == "1" | model_ver == "1.3") %>%
  filter(first_selection == "free") %>%
  # mutate(percent = 100*(fixed-free)/(free)) %>%
  ggplot(aes(dataset, abs_err, fill = as.factor(model_ver))) + 
  geom_col(position = position_dodge(), alpha = 0.75) +
  theme_bw() +
  #facet_grid(~first_selection, scales = "free_y") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(fill="model version") +
  ylab("absolute error") -> rl_z_a

rl_z %>%
  filter(model_ver == "1.3") %>%
  group_by(model_ver, first_selection) %>%
  summarise(abs_err = mean(abs_err)) %>%
  ggplot(aes(first_selection, abs_err)) + 
  geom_col(position = position_dodge(), alpha = 0.75) +
  theme_bw() +
  #facet_grid(~first_selection, scales = "free_y") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("absolute error") +
  xlab("first selection") -> rl_z_b

rl_z2 <- rl %>%
  ungroup() %>%
  mutate(abs_err = abs(observed - predicted), .keep = "unused") %>%
  #reframe(scale_abs_err = scale(abs_err, center = FALSE), .by = c(dataset, statistic)) %>%
  mutate(model_ver = rl$model_ver) %>%
  mutate(first_selection = if_else(str_detect(model_ver, "f"), "fixed", "free"),
         model_ver = str_remove(model_ver, "v|f"),
         model_ver = str_replace(model_ver, "_", "."),
         model_ver = as.numeric(model_ver)) %>%
  group_by(statistic, model_ver, first_selection) %>%
  summarise(abs_err = mean(abs_err)) 
           
rl_z2 %>%
  filter(model_ver == "1.3") %>%
  group_by(statistic, first_selection) %>%
  summarise(abs_err = mean(abs_err)) %>%
  ggplot(aes(statistic, abs_err, fill = as.factor(first_selection))) + 
  geom_col(position = position_dodge(), alpha = 0.75) +
  theme_bw() +
  #facet_grid(~first_selection, scales = "free_y") +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  ylab("absolute error") +
  labs(fill="first selection")  -> rl_z_c


rl_z_b + rl_z_c / rl_z_a

ggsave("figs/fig10_fixed_v_free.pdf", width = 8, height = 5)

