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
  
  ds %>% mutate(item_class = factor(item_class)) -> ds
  
  # plot basic trial
  plt <- ggplot(data = ds, aes(x, y)) + 
    geom_path(data = df, colour = "grey", group = 1) +
    geom_point(size = 5, aes(colour = item_class, shape = item_class)) +
    geom_label_repel(data = df, aes(label = found), size = 3) + 
    scale_colour_manual(values = c(18, 15, 3, 4)) + 
    scale_shape_manual(values = c(15, 19, 3, 4))
  
  # if segLabel has been left unspecified, add a simple geom_path to connect selected items
  if (is_null(segLabel)) {
    
    plt <- plt 
    
  } else { # use geom_labelsegment..
    
    #first we have to rearrange things a little
    df %>% mutate(x0 = lag(x), y0 = lag(y)) -> df
    
    # turn labels into characeters and repace NA with a dash (-)
    df[[segLabel]] = as.character(round(df[[segLabel]], 2))
    df[[segLabel]][is.na(df[[segLabel]])]  = "-"
    
    plt <- plt + geom_textsegment(data = df, 
                                  aes(x = x0, y = y0, 
                                      xend = x, yend = y,
                                      label = .data[[segLabel]]), 
                                  colour = "grey")
  }
  
  plt <- plt + coord_equal() + 
    theme(axis.title = element_blank(),
          axis.ticks  = element_blank(),
          axis.text = element_blank(),
          legend.position = "none",
          plot.background = element_rect(fill='darkgrey', colour='black'),
          panel.grid = element_blank())
  
  if (!is.na(filename)) {
    ggsave(filename)
  }
  
  return(plt)
}


dataset <- "clarke2022qjep"

# read in data
d <- import_data(dataset)

plt1 <- plot_a_trial(d$stim, d$found, trial = 41)

plt2 <- plot_a_trial(d$stim, d$found, trial = 1263)

plt3 <- plot_a_trial(d$stim, d$found, trial = 201)

plt1 + plt2 + plt3
