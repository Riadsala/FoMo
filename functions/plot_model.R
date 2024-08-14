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


plot_model_human_iisv_comparison <- function(pred, df, iisv_emp = NULL) {
  
  # Create a plot to compare our model to human/training data
  # in terms of  inter-item selection amplitude and direction
  
  iisv_sim <- get_iisv_over_trials(pred$sim)
  
  if (is.null(iisv_emp)) iisv_emp <- get_iisv_over_trials(df$found) 
  
  if (!("person" %in% names(iisv_emp))) iisv_emp$person = 1
  
  #################################################################
  # create distance plot 
  iisv_emp %>% 
    filter(found > 1) %>%
    group_by(found, condition) %>%
    summarise(distance = median(sqrt(d2)), .groups = "drop") -> emp 
  
  iisv_sim %>%
    filter(found > 1) %>%
    group_by(found, condition, .draw) %>%
    summarise(distance = median(sqrt(d2)), .groups = "drop_last") %>%
    median_hdci(distance, .width = c(0.53, 0.97)) -> sim
  
  ggplot(emp, aes(found, distance)) + 
    geom_path(data = emp, 
              aes(colour = condition)) +
    geom_ribbon(data = sim, 
                aes(ymin = .lower, ymax = .upper, 
                    fill = condition,
                    group = interaction(.width, condition)), alpha = 0.33) -> plt_amp
  
  if (length(unique(iisv_sim$condition)) == 1) {
    
    plt_amp <- plt_amp + theme(legend.position = "none")
  }
  
  #################################################################
  # create direction plot
  
  # first, we need to pad our data to take 0/2pi into account
  bind_rows(iisv_emp, 
            iisv_emp %>% mutate(theta = theta + 2*pi)) %>%
    filter(found > 1) -> iisv_emp2
  
  bind_rows(iisv_sim, 
            iisv_sim %>% mutate(theta = theta + 2*pi)) %>%
    filter(found > 1) -> iisv_sim2
  
  
  pi_labels <- c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi))
  
  iisv_emp2 %>%
    ggplot(aes(theta, group = person)) + 
    geom_line(stat = "density", bw = 0.1) +
    geom_line(data = iisv_sim2, aes(group = .draw), 
              alpha = 0.25, stat = "density", bw = 0.1, colour = "red") +
    coord_cartesian(xlim = c(0, 2*pi)) +
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    scale_x_continuous("inter-item directions", breaks = seq(0, 2*pi, pi/2), 
                       labels = pi_labels) -> plt_wave
  
  rm(iisv_sim2, iisv_emp2)
  
  #################################################################
  # create rel direction plot
  
  # first, mirror to fix density plot around 0 
  bind_rows(iisv_emp, 
            iisv_emp %>% mutate(psi = -psi)) %>%
    filter(found > 2) -> iisv_emp2
  
  bind_rows(iisv_sim, 
            iisv_sim %>% mutate(psi = -psi)) %>%
    filter(found > 2) -> iisv_sim2
  
  iisv_emp2 %>%
    ggplot(aes(psi, group = person)) + 
    geom_line(stat = "density", bw = 0.1) +
    geom_line(data = iisv_sim2, aes(group = .draw), 
              alpha = 0.25, stat = "density", bw = 0.1, colour = "red") + 
    coord_cartesian(xlim = c(0, 1)) -> plt_psi
  
  rm(iisv_sim2, iisv_emp2)
  
  
  #################################################################
  # output plots
  
  plt <- plt_amp + plt_wave + plt_psi
  
  return(plt)
}





plot_model_accuracy <- function(pred) {
  
  n_targets <- max((pred$acc$found))
  
  baseline <- tibble(found = 1:n_targets, accuracy = 1/((n_targets + 1) - found))
  
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