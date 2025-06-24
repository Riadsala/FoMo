library(tidyverse)
library(patchwork)
library(circular)

source("../functions/import_data.R")
source("../functions/prep_data.R")
# source("../functions/compute_summary_stats.R")
source("../functions/plot_model.R")
source("../functions/plot_data.R")
source("../functions/post_functions.R")
source("../functions/sim_foraging_data.R")

options(mc.cores = 1, digits = 2)
# set global ggplot theme
theme_set(theme_bw())
sf <- "../examples/1_fit_models/scratch"

# select dataset
ds <- "clarke2022qjep"

# customised plot_a_trial function
plot_a_trial <- function(ds, df, 
                         trial = NA, segLabel = NULL,
                         filename = NA,
                         draws_to_plot = 1,
                         pred) {
  
  # This function plots a trial: points indicate all items and a path 
  # joins the items up in the order in which they where selected.
  # Path vertices are labelled to indicate the order in which they were selected
  # Path segments may be labelled with some feature (a column in df).
  
  trl = trial
  ds <- filter(ds, trial == trl)
  df <- filter(df, trial == trl)
  
  # split stimuli into two tibbles so we can plot each with different colours
  ds1 <- filter(ds, item_class == 1)
  ds2 <- filter(ds, item_class == 2)
  
  ds %>% mutate(item_class = factor(item_class)) -> ds
  
  # add predictions
  pred %>% filter(trial_p == trl, .draw %in% c(1,2,3)) -> predf
  
  my_theme <- theme(legend.position = "none", 
                    axis.title = element_blank(),
                    axis.ticks  = element_blank(),
                    axis.text = element_blank(),
                    # plot.background = element_rect(fill='darkgrey', colour='black'),
                    panel.grid = element_blank(),
                    panel.background =  element_rect(fill='darkgrey', colour='darkgrey'))

  # plot human path
  plt_human <- ggplot(ds, aes(x, y)) + 
    geom_path(data = df)  +
    geom_point(data = ds1, size = 5, aes(shape = factor(item_class)), color = "#A1CAF1") +
    geom_point(data = ds2, size = 5, aes(shape = factor(item_class)), color = "#BE0032") +
    geom_text(data = df, aes(label = found), size = 2.5) + 
    scale_shape_manual(values = c(19, 19, 3, 4))+ 
    coord_equal() + 
    scale_linewidth(guide = "none") + 
    # scale_shape(guide = "none") +
    my_theme
  
  
  # plot model simulations
  plt_model <- ggplot(data = ds, aes(x, y)) + geom_path(data = predf, 
                  aes(x, y),
                  alpha = 0.5) + 
    geom_point(data = ds1, size = 5, aes(shape = factor(item_class)), color = "#A1CAF1") +
    geom_point(data = ds2, size = 5, aes(shape = factor(item_class)), color = "#BE0032") +
    geom_text(data = predf, aes(label = found), size = 2.5) + 
    scale_color_paletteer_d("wesanderson::Chevalier1") +
    scale_shape_manual(values = c(19, 19, 3, 4)) + coord_equal() + 
    facet_grid(first_selection ~ .draw) +
     scale_linewidth(guide = "none") + 
    # scale_shape(guide = "none") +
    my_theme
  
  pltout <- plt_human + plt_model + plot_layout(widths = c(2,3))

  return(pltout)
}

# read in data and predictions
d <- import_data(ds)
folder <- paste0(sf, "/post/", ds, "/")

pred13 <- readRDS(paste0(folder, "pred_1_4.rds"))

bind_rows(
  pred13$trialwise %>% mutate(first_selection = "unconstrained"),
  pred13$trialwise_firstfixed %>% mutate(first_selection = "constrained")) -> pred

plot_a_trial(d$stim, d$found, trial = 17, pred = pred)

ggsave("figs/fig8_examples.pdf", width = 15, height = 6)

