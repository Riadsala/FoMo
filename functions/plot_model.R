library(tidyverse)
library(patchwork)
library(tidybayes)
library(hrbrthemes)

# plotting our foraging model

# plot_model_fixed() plots the fixed effects (group averages) 
# plot_post_prior() - why is this called this? 
# plot_traceplots(m) plots the traceplots for the core params in our model 
# plot_beta_comp() details needed
# plot_init_sel() - plots the initial selection component of the model

# plot_init_sel() needs a tidy up?
# MORE sample_beta() to post_functions ???
# plot_model_spatial () will need fixed at some point

theme_set(theme_ipsum())
options(ggplot2.discrete.colour = ggthemes::ptol_pal()(2),
        ggplot2.discrete.fill = ggthemes::ptol_pal()(2))


plot_model_human_rl_comparison <- function(rl) {
  
  #################################################################
  # plot scatter plots for mean max run length and mean num runs
  
  rl %>% 
    pivot_longer(c("max_run_length", "n_runs"), names_to = "statistic") %>%
    pivot_wider(names_from = "data") %>%
    ggplot(aes(observed, simulated, colour = condition)) + 
    geom_point() +
    facet_wrap(~statistic) -> plt
  
  return(plt)
  
}

plot_model_human_iisv_comparison <- function(iisv) {
  
  # Create a plot to compare our model to human/training data
  # in terms of  inter-item selection amplitude and direction
  
  #################################################################
  # create distance plot
  ggplot(iisv, aes(d2, colour = data, group = interaction(.draw, data))) + 
    geom_density(alpha = 0.33) +
    scale_colour_viridis_d() +
    facet_wrap(~condition) -> plt_amp
  
  #################################################################
  # create direction plot
  pi_labels <- c("0", expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi))
  
  iisv %>%
    ggplot(aes(theta, colour = data, group = interaction(.draw, data))) + 
    geom_line(stat = "density", bw = 0.1, alpha = 0.33,
              linewidth = 2) +
    coord_cartesian(xlim = c(-pi, pi)) +
    theme(legend.position = "none") +
    scale_colour_viridis_d() +
    scale_x_continuous("inter-item directions", breaks = seq(0, 2*pi, pi/2), 
                       labels = pi_labels) -> plt_wave
  
  #################################################################
  # create rel direction plot
  
  # first, mirror to fix density plot around 0 
  bind_rows(iisv, 
            iisv %>% mutate(psi = -psi)) %>%
    filter(found > 2) -> iisv
  
  iisv %>%
    ggplot(aes(psi,  colour = data, group = interaction(.draw, data))) + 
    geom_line(stat = "density", bw = 0.1,
               linewidth = 2, alpha = 0.5)  + 
    scale_colour_viridis_d() +
    coord_cartesian(xlim = c(0, 1)) -> plt_psi
  
  
  #################################################################
  # output plots

  plt <- plt_amp + plt_psi + plt_wave + plot_layout(guides = "collect")

  return(plt)
}

plot_model_accuracy <- function(acc) {
  
  # acc is a dataframe of .draws, as computed by  summarise_acc() in post_functions.R
  
  #calculate chance baseline
  n_targets <- max((acc$found))
  baseline <- tibble(found = 1:n_targets, accuracy = 1/((n_targets + 1) - found))
  
  # work out grouping factors
  groups <- names(acc)[!(names(acc) %in% c(".draw", "accuracy"))]
  
  # calculate HPDI
  acc %>%
    group_by_at(groups) %>%
    median_hdci(accuracy) -> acc_hpdi
  
   acc_hpdi %>%
    ggplot(aes(found, accuracy)) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = condition), 
                alpha = 0.5) + 
    geom_line(aes(colour = condition)) + 
    geom_path(data = baseline, linetype = 2) -> plt
  
  if ("split" %in% names(acc)) {
    plt <- plt + facet_wrap(~split)
  } 
  
  return(plt)
}

plot_model_fixed <- function(post, gt=NULL, clist=NULL)
{
  
  # sort out groudtruth params if passed in as sim params
  if ("foraging" %in% names(gt)) {
    
    gt <- gt$foraging
    
    gt <- list(b_a = qlogis(filter(gt, param == "bA")$mu[[1]][1]),
               b_stick = filter(gt, param == "bS")$mu[[1]],
               rho_delta = filter(gt, param == "rho_delta")$mu[[1]],
               rho_psi = filter(gt, param == "rho_psi")$mu[[1]])
    
    print(gt)
    
  }
  
  my_widths <- c(0.53, 0.97)
  
  # create a plot for each parameter
  plts <- map(post$params, plt_post_prior, 
              post = post$fixed, prior = post$prior, 
              gt = gt, clist = clist)
  
  # assemble the plots!
  plt <- wrap_plots(plts, nrow = 1) + 
    plot_layout(guides = "collect", axis_titles = "collect")
  
  return(plt)
}

plot_cts_params <- function(post,  gt=NULL, clist=NULL)
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
    as_tibble() %>%
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
  
  plts <- map(params, create_weight_plot, post = post)
  
  plt <- wrap_plots(plts) + plot_layout(guides = "collect")
  
  return(plt)
  
}

create_weight_plot <- function(param, x1, x2, post) {
  
  x1 <- 0
  
  if (param == "rho_psi") {
    x2 <- 1  
  } else {
    x2 <- 2
  }
  
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