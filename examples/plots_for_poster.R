library(tidyverse)
library(patchwork)
library(tidybayes)
library(hrbrthemes)
library(paletteer)

ds <- "bhat2025"
scratch_folder <- "1_fit_models/scratch"

### accuracy plot

files <- dir(paste0(scratch_folder, "/post/", ds))
files <- files[str_detect(files,"pred")]

d <- tibble()

for (ff in files) {
  
  a <- readRDS(paste0(scratch_folder, "/post/", ds, "/", ff))
  
  a$itemwise$model_ver <- a[[4]]
  
  d <- bind_rows(d, a$itemwise)
  
}

# for now, only want 1.0 and 1.3

d <- filter(d, model_ver == "1_0" | model_ver == "1_3")

facet_names <- c(
  conj_grid = "conjunction cardinal",
  `conj_rot` = "conjunction rotated",
  `feat_grid` = "feature cardinal",
  `feat_rot` = "feature rotated"
)

xlabs <- c("original", "new")

d %>% 
  group_by(model_ver, condition, .draw) %>%
  summarise(accuracy = mean(model_correct),
            .groups = "drop") %>%
  ggplot(aes(model_ver, accuracy)) +
  stat_interval(alpha = 0.5) + 
  facet_wrap(~condition, labeller = as_labeller(facet_names)) +
  xlab("model version") +
  scale_x_discrete(labels= xlabs) +
  theme_bw(base_size = 18) +
  theme(legend.position = "none",  
        panel.grid.major = element_blank()) 
  

ggsave("plots/accuracy_plot.png", width = 8, height = 5)

### accuracy by person plot

v1 <- "1_0"
v2 <- "1_3"

acc <- tibble()

folder <- paste0(scratch_folder, "/post/", ds, "/")
  
acc1 <- readRDS(paste0(folder, "pred_", v1, ".rds"))$itemwise %>%
  mutate(version = v1) %>%
  filter(found > 1, found < 40) %>%
  group_by(version, person, condition, .draw) %>%
  summarise(accuracy = mean(model_correct),
            .groups = "drop_last") %>%
  median_hdci(accuracy)
  
acc2 <- readRDS(paste0(folder, "pred_", v2, ".rds"))$itemwise %>%
  mutate(version = v2) %>%
  filter(found > 1, found < 40) %>%
  group_by(version, person, condition, .draw) %>%
  summarise(accuracy = mean(model_correct),
            .groups = "drop_last") %>%
  median_hdci(accuracy)
  
acc <- bind_rows(acc, 
                   bind_rows(acc1, acc2) %>% mutate(dataset = ds))

rm(acc1, acc2)

acc %>% 
  unite(accuracy, .lower, accuracy, .upper) %>%
  pivot_wider(names_from = "version", values_from = "accuracy") %>%
  separate(v1, c("xmin", "x", "xmax"), "_", convert = TRUE) %>%
  separate(v2, c("ymin", "y", "ymax"), "_", convert = TRUE) %>%
  mutate(improvement = cut(y - x, breaks = seq(-.025, .50, 0.05))) %>%
  separate(improvement, into = c("l", "u"), sep = ",") %>%
  mutate(improvement = factor((parse_number(l)+parse_number(u))/2),
         .keep = "unused")-> acc

acc %>% drop_na() %>% ggplot(aes(x, y, 
                   xmin = xmin, xmax = xmax,
                   ymin = ymin, ymax = ymax,
                   colour = improvement)) +
  geom_point(alpha = 0.9) + 
  geom_errorbar(alpha = 0.25) +
  geom_errorbarh(alpha = 0.25) + 
  geom_abline(linetype = 2) + 
  coord_equal() +
  scale_color_paletteer_d("MoMAColors::Abbott", direction = -1) +
  facet_wrap(~ condition, labeller = as_labeller(facet_names)) +
  theme_bw(base_size = 26) +
  theme(legend.position = "none", panel.grid.major = element_blank()) + 
  xlab("original FoMo model accuracy") +
  ylab("new FoMo model accuracy") 
  
ggsave("plots/accuracy_by_person_plot.png", width = 8, height = 5)

### looking at a case where model has improved things

acc <- tibble()

acc1 <- readRDS(paste0(folder, "pred_", v1, ".rds"))$itemwise %>%
  mutate(version = v1) %>%
  filter(found > 1, found < 40) %>%
  group_by(version, person, condition, .draw) %>%
  summarise(accuracy = mean(model_correct),
            .groups = "drop_last") %>%
  median_hdci(accuracy)

acc2 <- readRDS(paste0(folder, "pred_", v2, ".rds"))$itemwise %>%
  mutate(version = v2) %>%
  filter(found > 1, found < 40) %>%
  group_by(version, person, condition, .draw) %>%
  summarise(accuracy = mean(model_correct),
            .groups = "drop_last") %>%
  median_hdci(accuracy)

acc <- bind_rows(acc, 
                 bind_rows(acc1, acc2) %>% mutate(dataset = ds))

rm(acc1, acc2)

acc_of_interest <- acc %>%
  filter(condition == "feat_grid")

# person 1 and person 13 seem to show quite a lot of improvement (used trial 303)
# model is generally not very good for person 27. 18 also doesn't improve much (used trial 420)

source("../functions/import_data.R")

d <- import_data("bhat2025")

ds <- d$stim
df <- d$found

d_for_plot <- filter(df, condition == "feat_grid", person == "18")

#plot_a_trial(ds, df, trial = 303)

trial = 420

trl = trial
ds <- filter(ds, trial == trl)
df <- filter(df, trial == trl)

ds %>% mutate(item_class = factor(item_class)) -> ds

plt <- ggplot(data = ds, aes(x, y)) + 
  geom_point(size = 8, aes(colour = item_class, shape = item_class)) +
  ggrepel::geom_label_repel(data = df, aes(label = found), size = 8) + 
  scale_colour_manual(values = c(18, 15, 3, 4)) + 
  scale_shape_manual(values = c(15, 19, 3, 4))
  
plt <- plt + geom_path(data = df, colour = "grey", linewidth = 2, group = 1) 
  
plt <- plt + coord_equal() + 
  theme(axis.title = element_blank(),
        axis.ticks  = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.background = element_rect(fill='darkgrey', colour='black'),
        panel.grid = element_blank())

plt
