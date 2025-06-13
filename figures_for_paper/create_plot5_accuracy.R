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

datasets <- c("hughes2024rsos", "clarke2022qjep")   


v1 <- "1_0"
v2 <- "1_3"

plt1 <- plot_model_accuracy_comparison(c("clarke2022qjep", "hughes2024rsos"), v1, v2, scratch_folder = sf)

 acc <- tibble()

for (ds in datasets) {
  folder <- paste0(sf, "/post/", ds, "/")
  pred <- readRDS(paste0(folder, "pred_1_3.rds"))
  
  pred$itemwise %>% 
    group_by(condition, found, .draw, person) %>%
    summarise(person_acc = mean(model_correct), .groups = "drop_last") %>%
    summarise(accuracy = mean(person_acc), .groups = "drop_last") %>% 
    mutate(dataset = ds) %>%
    bind_rows(acc) -> acc
}

acc %>% group_by(dataset) %>% summarise(n_targets = max(found)+1)
 
baseline <- bind_rows(
   tibble(found = 1:40,
          accuracy = 1/((40 + 1) - found),
          dataset = "clarke2022qjep"),
   tibble(found = 1:20,
          accuracy = 1/((20 + 1) - found),
          dataset = "hughes2024rsos")
 )
  

ggplot(acc, aes(found, accuracy)) +
  stat_lineribbon(aes(fill = condition), alpha = 0.5, linewidth = 0) +
  geom_path(data = baseline) + 
  facet_wrap(~dataset, scales = "free", nrow = 2) -> plt2

plt1 + plt2 + plot_layout(widths = c(2,1))

ggsave("figs/fig5_acc.pdf", width = 10, height = 4)  
