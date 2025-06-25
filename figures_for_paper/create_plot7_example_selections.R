library(tidyverse)
library(patchwork)
library(circular)

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
sf <- "../examples/1_fit_models/scratch"

# select dataset
ds <- "clarke2022qjep"

# customised plot_a_trial function
plot_a_trial <- function(ds, df, 
                         trial = NA, segLabel = NULL,
                         filename = NA,
                         draws_to_plot = 1,
                         pred, ff) {
  
  # This function plots a trial: points indicate all items and a path 
  # joins the items up in the order in which they where selected.
  # Path vertices are labelled to indicate the order in which they were selected
  # Path segments may be labelled with some feature (a column in df).
  
  trl = trial
  ds <- filter(ds, trial %in% trl)
  df <- filter(df, trial %in% trl)
  
  ds %>% mutate(item_class = factor(item_class)) -> ds
  
  # plot basic trial
  plt <- ggplot(data = ds, aes(x, y)) + 
    geom_path(data = df, colour = "grey80", group = 1) -> plt
  
  df %>% filter(found == ff-1) %>%
    rename(x0 = x, y0 = y) %>%
    select(trial, x0, y0) -> dff
  
  # add predictions
  pred %>% filter(found == ff, trial_p %in% trl ) %>%
    group_by(P, model, trial_p) %>%
    summarise(n = n()) %>%
    rename(id = "P", trial = "trial_p") %>%
    left_join(ds, by = join_by(id, trial)) %>%
    left_join(dff) %>%
    filter(n > 2) -> predf
  
  plt + geom_segment(data = predf, 
                     aes(x=x0, y=y0, xend = x, yend = y, linewidth = n),
                     alpha = 0.5, colour = "darkviolet")  -> plt
  
  print(predf %>% group_by(model) %>% summarise(prop_acc_for = sum(n)/max(pred$.draw)))
  
  ds1 <- filter(ds, item_class == 1)
  ds2 <- filter(ds, item_class == 2)
  
   plt + geom_point(data = ds1, size = 5, aes(shape = item_class), color = "#A1CAF1") +
     geom_point(data = ds2, size = 5, aes(shape = item_class), color = "#BE0032") +
    geom_text(data = df, aes(label = found), size = 2.5) + 
     scale_color_paletteer_d("wesanderson::Chevalier1") +
    scale_shape_manual(values = c(19, 19, 3, 4)) -> plt
  
  plt <- plt + coord_equal() + 
    facet_grid(trial~model) +
    scale_linewidth(guide = "none") + 
    scale_shape(guide = "none") +
    theme(axis.title = element_blank(),
          axis.ticks  = element_blank(),
          axis.text = element_blank(),
          # plot.background = element_rect(fill='darkgrey', colour='black'),
          panel.grid = element_blank(),
          panel.background =  element_rect(fill='darkgrey', colour='darkgrey'))

  return(plt)
}

# read in data and predictions
d <- import_data(ds)
folder <- paste0(sf, "/post/", ds, "/")

pred10 <- readRDS(paste0(folder, "pred_1_0.rds"))$itemwise %>% 
  mutate(model = "1.0")
pred13 <- readRDS(paste0(folder, "pred_1_3.rds"))$itemwise %>% 
  mutate(model = "1.3")

plot_a_trial(d$stim, d$found, trial = c(53, 1273, 211), pred = bind_rows(pred10, pred13), ff = 8)


