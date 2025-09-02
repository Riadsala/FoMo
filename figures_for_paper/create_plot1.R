library(tidyverse)
library(patchwork)

source("../functions/import_data.R")
source("../functions/compute_summary_stats.R")

options(mc.cores = 1, digits = 2)

# set global ggplot theme
theme_set(theme_bw())

plot_a_trial <- function(ds, df, 
                         trial = NA, segLabel = NULL,
                         filename = NA,
                         draws_to_plot = 1) {
  
  # This function plots a trial: points indicate all items and a path 
  # joins the items up in the order in which they where selected.
  # Path vertices are labelled to indicate the order in which they were selected
  # Path segments may be labelled with some feature (a column in df).
  
  trl = trial
  ds <- filter(ds, trial == trl)
  df <- filter(df, trial == trl)
  
  ds %>% mutate(item_class = factor(item_class)) -> ds
  
  # plot basic trial
  plt <- ggplot(data = ds, aes(x, y)) + 
    geom_path(data = df, colour = "grey80", group = 1) +
    geom_point(size = 5, aes(colour = item_class, shape = item_class)) +
    geom_text(data = df, aes(label = found), size = 2.5) + 
    scale_colour_manual(values = c("#F0F8FF", "#BE0032", 3, 4)) + 
    scale_shape_manual(values = c(19, 19, 3, 4))
  
  plt <- plt + coord_equal() + 
    theme(axis.title = element_blank(),
          axis.ticks  = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          # plot.background = element_rect(fill='darkgrey', colour='black'),
          panel.grid = element_blank(),
          panel.background =  element_rect(fill='darkgrey', colour='darkgrey'))
  
  return(plt)
}


dataset <- "clarke2022qjep"

# read in data
d <- import_data(dataset)

plt1 <- plot_a_trial(d$stim, d$found, trial = 41)

plt2 <- plot_a_trial(d$stim, d$found, trial = 1263)

plt3 <- plot_a_trial(d$stim, d$found, trial = 201)

plt1 + plt2 + plt3

ggsave("figs/fig1_example_trials.pdf", width = 12, height = 3.2)


