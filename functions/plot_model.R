library(tidyverse)
library(patchwork)
library(tidybayes)

# plotting our foraging model

# plot_model_fixed() plots the fixed effects (group averages) 
# plot_post_prior() - why is this called this? 
# plot_traceplots(m) plots the traceplots for the core params in our model 
# plot_beta_comp() details needed
# plot_init_sel() - plots the initial selection component of the model

# plot_init_sel() needs a tidy up?
# MORE sample_beta() to post_functions ???
# plot_model_spatial () will need fixed at some point

theme_set(theme_minimal())
options(ggplot2.discrete.colour = ggthemes::ptol_pal()(2),
        ggplot2.discrete.fill = ggthemes::ptol_pal()(2))


plot_model_human_comparison <- function(pred, df) {
  
  # Create a plot to compare our model to human/training data
  # in terms of run statistics, inter-target length etc.
  
  d <- get_inter_sel_info_over_trials(post$sim)
  
  e <- get_inter_sel_info_over_trials(df$found) %>%
    group_by(person, found) %>% 
    summarise(empirical = mean(d2),
              .groups = "drop")
  
  d %>% group_by(person, found) %>% 
    summarise(simulated = mean(d2)) %>%
    full_join(e, by = join_by(person, found)) %>%
    filter(found != 1) %>%
    pivot_longer(c(simulated, empirical), values_to = "d2") %>%
    ggplot(aes(found, sqrt(d2), colour = name, group = interaction(name, person))) + 
    geom_path(alpha = 0.5) +
    geom_smooth(aes(group = name), se=FALSE, linetype = 2) + 
    scale_x_continuous("nth item found") + 
    scale_y_continuous("inter-target distance") +
    theme(legend.title = element_blank()) -> plt_amp
  
  ############## now do run statistics
  
  d <- get_run_info_over_trials(post$sim)
  
  e <- get_run_info_over_trials(df$found) %>%
    group_by(person, condition) %>% 
    summarise(n_runs = median(n_runs),
              max_run_length = median(max_run_length),
              .groups = "drop") %>%
    mutate(type = "empirical") 
  
  d %>% group_by(person, condition) %>% 
    summarise(n_runs = median(n_runs),
              max_run_length = median(max_run_length),
              .groups = "drop") %>%
    mutate(type = "simulated") %>%
    full_join(e, by = join_by(person, n_runs, max_run_length, type, condition)) %>%
    pivot_longer(c(max_run_length, n_runs), names_to = "stat") %>%
    pivot_wider(names_from = "type", values_from ="value") %>%
    ggplot(aes(empirical, simulated, 
               colour = condition, shape = condition)) + 
    geom_abline(linetype = 2) + 
    geom_point(size = 2) +
    facet_wrap(~stat, scales = "free") -> plt_runs
  
  
  plt <- plt_amp / plt_runs
  
  return(plt)
}





plot_model_accuracy <- function(pred) {
  
  n_targets <- max((pred$acc$found))
  
  baseline <- tibble(found = 1:n_targets, accuracy = 1/(41 - found))
  
  pred$acc %>% group_by(found, .draw) %>%
    summarise(accuracy = mean(model_correct, .groups ="last")) %>%
    median_hdci(accuracy) %>%
    ggplot(aes(found, accuracy)) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.5) + 
    geom_path() + 
    geom_path(data = baseline, linetype = 2)
}

plot_model_fixed <- function(post, gt=NULL, clist=NULL)
{
  
  my_widths <- c(0.53, 0.97)
  
  # create a plot for each parameter
  plts <- map(post$params, plt_post_prior, 
              post = post$fixed, prior = post$prior, 
              gt = gt, clist = clist)
  
  # assemble the plots!
  plt <- wrap_plots(plts, nrow = 1) + 
    plot_layout(guides = "collect")
  
  return(plt)
}

plt_post_prior <- function(post, prior, var, gt=NULL, clist=NULL) {
  
  # function to plot the posterior against the prior. 
  # gt allows us to mark up the groundtruth (if available)
  # clist allows us to specify a list of conditions to use in different ways
  
  if (is.null(clist)) {
    fill_cond <- "condition"
  } else {
    fill_cond <- clist$fill
  }
  
  prior_var = paste("prior", var, sep = "_")
  
  # get prior HDPI
  prior %>% 
    median_hdci(get(prior_var), .width = c(0.53, 0.97)) -> prior_hpdi
  
  post %>% 
    ggplot() + 
    geom_rect(data = prior_hpdi,
              aes(ymin = -Inf, ymax = Inf, xmin = .lower, xmax = .upper), 
              fill = "grey", alpha = 0.25) +  
    geom_density(aes(get(var), fill = !!sym(fill_cond)), alpha = 0.5) +
    scale_x_continuous(var) -> plt
  
  if (!is.null(gt)) {
    
    plt <- plt + geom_vline(xintercept = gt[[var]], linetype = 2, colour = "darkred")
    
  }
  
  if (!is.null(clist)) {
    plt <- plt + facet_wrap(as.formula(paste("~", clist$facet_cond)))
  }
  
  return(plt)
  
}

plot_model_random <- function(post) 
{
  
  post$random %>%
    pivot_longer(starts_with("u"), names_to = "param") %>%
    group_by(person, condition, param) %>%
    median_hdci(value) -> d_hpdi
  
  d_hpdi %>% ggplot(aes(person, ymin = .lower, y = value, ymax = .upper, color = condition)) +
    geom_hline(linetype = 2, yintercept = 0) + 
    geom_linerange() +
    facet_wrap(~param, scales = "free") -> plt
  
  return(plt)
  
}

plot_model_weights <- function(post, params) {
  
  plts <- pmap(params, create_weight_plot, post = post)
  
  plt <- wrap_plots(plts) + plot_layout(guides = "collect")
  
}

create_weight_plot <- function(param, x1, x2, post) {
  
  x <- seq(x1, x2, (x2-x1)/100)
  
  rho <- post$fixed %>% select(condition, "rho" = {param}) 
  dout <- pmap_df(rho, neg_exp, x=x)
  
  plt <- ggplot(dout, aes(x = x, y= p, fill = condition)) + 
    stat_lineribbon(alpha = 0.4) +
    scale_x_continuous(param)
  
  return(plt)
  
} 

neg_exp <- function(condition, rho, x) {
  
  return(tibble(condition = condition, 
                rho = rho, 
                x = x, 
                p = exp(-rho*x)))
}

plot_traceplots <- function(m) 
{
  
  bayesplot::mcmc_trace(m$draws(), pars = c("bA", "b_stick"))
  
}