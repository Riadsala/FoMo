d <- readRDS("scratch/d_1_0.rds")
fit <- readRDS("scratch/unknown_all_1_0.model")

post <- extract_post(fit, d, multi_level = TRUE)

plot_model_fixed(post, gt = list(b_a = qlogis(item_class_weights[[1]][1]),
                                 b_stick = b_stick,
                                 rho_delta = rho_delta,
                                 rho_psi = rho_psi))

pred <- summarise_postpred(fit, d, 
                             get_sim = TRUE, draw_sample_frac=0.002)


acc <- compute_acc(pred$acc) %>% mutate(model = "FoMo 1.0")

# compute empirical run statistics
rle <- get_run_info_over_trials(d$found) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# compute simulated run statistics
rlp <- get_run_info_over_trials(pred$sim) %>%
  group_by(person, condition) %>%
  summarise(max_run_length = mean(max_run_length))

# bind everything together
bind_rows(rle %>% mutate(x = "observed"),
          rlp %>% mutate(x = "predicted")) %>%
  pivot_wider(names_from = "x", values_from = "max_run_length") -> rl

rl %>% ggplot(aes(observed, predicted)) + geom_point() +
  geom_abline()
