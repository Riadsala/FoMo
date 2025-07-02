library(tidyverse)
library(patchwork)

source("../functions/import_data.R")
source("../functions/compute_summary_stats.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

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

model_ver <- "v1_0"
dataset <- "hughes2024rsos"

# read in data
d <- import_data(dataset)

sf <- "../examples/1_fit_models/scratch"
folder <- paste0(sf, "/post/", dataset, "/")

#############################################################################
# compute & compare iisv statistics
#############################################################################

iisv <- read_csv(paste0("../examples/1_fit_models/scratch/post/", dataset, "/iisv_statistics.csv")) %>%
  rename(data = "z") %>%
  filter(data == "observed" | model_version == model_ver)

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

labels <- c(expression(-pi), expression(-pi/2), expression(0), expression(pi), expression(pi/2), expression(pi))

iisv %>% 
  filter(is.finite(theta)) %>%
  ggplot(aes(theta, fill = data)) + 
  geom_histogram(position = position_identity(),
                 breaks = seq(-pi, pi, pi/8), linewidth = 2,
                 alpha = 0.5) + 
  scale_x_continuous(expression(phi), 
                     breaks = c(-pi, -pi/2, 0, pi, pi/2, pi),
                     labels = labels) +
  paletteer::scale_fill_paletteer_d("lisa::BridgetRiley", direction = -1) +
  theme(legend.position = "bottom") -> plt_phi

ggsave("figs/figure2_absdir.pdf", width = 5, height = 4)


# compare 1.0 and 1.3

model_ver <- "f1_0"

iisv <- read_csv(paste0("../examples/1_fit_models/scratch/post/", dataset, "/iisv_statistics.csv")) %>%
  rename(data = "z") %>%
  filter(is.na(.draw) | .draw == 1) %>%
  filter(data == "observed" | model_version == model_ver)

iisv %>% 
  filter(is.finite(theta)) %>%
  ggplot(aes(theta, fill = data, colour = data)) + 
  geom_histogram(position = position_identity(),
                 breaks = seq(-pi, pi, pi/8), linewidth = 1,
                 alpha = 0.25) + 
  scale_x_continuous(expression(phi), 
                     breaks = c(-pi, -pi/2, 0, pi, pi/2, pi),
                     labels = labels) +
  paletteer::scale_fill_paletteer_d("lisa::BridgetRiley", direction = 1) +
  paletteer::scale_color_paletteer_d("lisa::BridgetRiley", direction = 1) +
  theme(legend.position = "bottom") -> plt_phi_v1

model_ver <- "f1_3"

iisv <- read_csv(paste0("../examples/1_fit_models/scratch/post/", dataset, "/iisv_statistics.csv")) %>%
  rename(data = "z") %>%
  filter(is.na(.draw) | .draw == 1)  %>%
  filter(data == "observed" | model_version == model_ver)

iisv  %>%
  filter(is.finite(theta)) %>%
  ggplot(aes(theta, fill = data, colour = data)) + 
  geom_histogram(position = position_identity(),
                 breaks = seq(-pi, pi, pi/8), linewidth = 1,
                 alpha = 0.25) + 
  scale_x_continuous(expression(phi), 
                     breaks = c(-pi, -pi/2, 0, pi, pi/2, pi),
                     labels = labels) +
  paletteer::scale_fill_paletteer_d("lisa::BridgetRiley", direction = 1) +
  paletteer::scale_color_paletteer_d("lisa::BridgetRiley", direction = 1) +
  theme(legend.position = "bottom") -> plt_phi_v3
 
plt_phi_v1 + plt_phi_v3 
