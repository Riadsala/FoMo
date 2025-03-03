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

model_ver <- "1_5"
dataset <- "clarke2022qjep"

# read in data
d <- import_data(dataset)

folder <- paste0("1_fit_models/scratch/post/", dataset, "/")


#############################################################################
# plot model comparison over models
#############################################################################
plot_models_accuracy(dataset)

# scatter plot of person acc by model
v1 <- "1_0"
v2 <- "1_3"

plot_model_accuracy_comparison <- function(dataset, v1, v2) {
  
acc1 <- readRDS(paste0(folder, "pred_train", v1, ".rds"))$itemwise %>%
  mutate(version = v1) %>%
  filter(found > 1, found < 40, split == "testing") %>%
  group_by(version, person, condition, .draw) %>%
  summarise(accuracy = mean(model_correct),
            .groups = "drop_last") %>%
  median_hdci(accuracy)

acc2 <- readRDS(paste0(folder, "pred_train", v2, ".rds"))$itemwise %>%
  mutate(version = v2) %>%
  filter(found > 1, found < 40, split == "testing") %>%
  group_by(version, person, condition, .draw) %>%
  summarise(accuracy = mean(model_correct),
            .groups = "drop_last") %>%
  median_hdci(accuracy)

acc <- bind_rows(acc1, acc2) 

rm(acc1, acc2)

acc %>% 
  unite(accuracy, .lower, accuracy, .upper) %>%
  pivot_wider(names_from = "version", values_from = "accuracy") %>%
  separate(v1, c("xmin", "x", "xmax"), "_", convert = TRUE) %>%
  separate(v2, c("ymin", "y", "ymax"), "_", convert = TRUE) %>%
  ggplot(aes(x, y, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_point(alpha = 0.75) + 
  geom_errorbar(alpha = 0.25) +
  geom_errorbarh(alpha = 0.25) + 
  geom_abline(linetype = 2) + 
  facet_wrap( ~ condition) + 
  coord_equal() +
  scale_x_continuous(paste0("FoMo v", str_replace(v1, "_", "."))) + 
  scale_y_continuous(paste0("FoMo v", str_replace(v2, "_", "."))) + 
  theme(panel.grid  = element_blank())


}

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
m <- readRDS(paste0("1_fit_models/scratch/models/", dataset, "/train", model_ver, ".model"))
post <- extract_post(m, d)
post_plt <- plot_model_fixed(post)

#############################################################################
# create plot
#############################################################################
acc_plt / post_plt

#############################################################################
# compute & compare run statistics
#############################################################################
pred <- readRDS(paste0(folder, "pred_train1_0.rds"))

rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs       = mean(n_runs),
            .groups = "drop")

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$trialwise %>% filter(.draw == 1)) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length),
            num_runs       = mean(n_runs))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) %>%
  pivot_longer(c(max_run_length, num_runs), names_to = "statistic") %>%
  pivot_wider(names_from = "x") -> rl

# plot
rl %>% ggplot(aes(predicted, observed, 
                  colour = condition, shape = condition)) + 
  geom_point() + 
  geom_abline(linetype = 2) + 
  facet_wrap(~statistic, scales = "free")

#############################################################################
# compute & compare iisv statistics
#############################################################################

rle <- get_iisv_over_trials(d$found) 

# compute simulated run statistics
rlp <- get_iisv_over_trials(pred$trialwise %>% filter(.draw == 1))

  # bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) -> rl

# plot
rl %>% ggplot(aes(predicted, observed, 
                  colour = condition, shape = condition)) + 
  geom_point() + 
  geom_abline(linetype = 2) + 
  facet_wrap(~statistic, scales = "free")
