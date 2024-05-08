library(tidybayes)
library(patchwork)

# plotting our foraging model

# plot_model_fixed() plots the fixed effects (group averages) 
# plot_post_prior() - why is this called this? 
# plot_traceplots(m) plots the traceplots for the core params in our model 
# plot_beta_comp() details needed
# plot_init_sel() - plots the inital selection component of the model

# plot_init_sel() needs a tidy up?
# MORE sample_beta() to post_functions ???
# plot_model_spatial () will need fixed at some point


options(ggplot2.discrete.colour = ggthemes::ptol_pal()(4),
        ggplot2.discrete.fill = ggthemes::ptol_pal()(4))

<<<<<<< HEAD
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

plot_model_fixed <- function(post, gt=NULL)
=======
plot_model_fixed <- function(post, m, d, cl = "A", gt=NULL, directions = FALSE)
>>>>>>> parent of a5f1521 (Joining in on the merge)
{
  # gt is a list with our groundtruth sim parameters. 
  post <- post$fixed
  
  m %>% 
    spread_draws(prior_bA, prior_b_stick, prior_rho_delta, prior_rho_psi) %>%
    ungroup() -> prior
  
  my_widths <- c(0.53, 0.97)
  
  # plot class weights
  post %>%
    ggplot() + 
    geom_rect(data = prior %>% 
                median_hdci(prior_bA, .width = c(0.53, 0.97)) %>%
                mutate(.lower = plogis(.lower), .upper = plogis(.upper)),
              aes(ymin = -Inf, ymax = Inf, xmin = .lower, xmax = .upper), 
              fill = "orange", alpha = 0.25) +
    geom_density(aes(plogis(bA), fill = condition), alpha = 0.5) +
    geom_vline(xintercept = 0.5, linetype = 2) +
    scale_x_continuous("class weights") +
    theme(legend.position = "bottom") -> plt_pA
  
  # plot stick-switch param
  post  %>%
    ggplot()  + 
    geom_rect(data = prior %>% 
                median_hdci(prior_b_stick, .width = c(0.53, 0.97)) %>%
                mutate(.lower = plogis(.lower), .upper = plogis(.upper)),
              aes(ymin = -Inf, ymax = Inf, xmin = .lower, xmax = .upper), 
              fill = "orange", alpha = 0.25) + 
    geom_density(aes(plogis(b_stick), fill = condition), alpha = 0.5) +
    geom_vline(xintercept = 0.5, linetype = 2) +
    scale_x_continuous("stick probability")  +
    theme(legend.position = "bottom") -> plt_pS
  
  # plot proximity and direction effects
  plt_prox    <- plt_post_prior(post, prior, "rho_delta", "proximity tuning", gt$rho_delta)
  
 # plt_mem    <- plt_post_prior(post, prior, "bM", "one back weighting")
  plt_rel_dir <- plt_post_prior(post, prior, "rho_psi", "direction tuning", gt$rho_psi)
 # plt_dir2 <- plt_post_prior(post, prior, "direction_bias", "Hori-Vert Pref") 
  
  
  if (!is.null(gt)) {
    # if we have groundtruth, annotate our plt
    plt_pA + geom_vline(xintercept = gt$pA, linetype = 1, colour = "red") -> plt_pA
    
    plt_pS + geom_vline(xintercept = plogis(gt$b_stick), linetype= 1, colour = "red") -> plt_pS
  }
  
  plt <- plt_pA + plt_pS + plt_prox + plt_rel_dir
  
  # add direction plot, if we're using that model
  if (directions) {
    m %>% recover_types(d) %>%
      spread_draws(theta[angle]) %>%
      mutate(angle = (angle-1)*90,
             angle = as_factor(angle)) %>%
      ungroup() %>%
      ggplot(aes(theta, fill = angle)) + 
      geom_density(alpha = 0.5) + 
      geom_vline(data =  tibble(angle = as_factor(0:1 * 90), theta = gt$theta), 
                 aes(xintercept = theta, colour = angle), linetype = 1) + 
      ggthemes::scale_color_ptol() -> plt_abs_dir_theta
    
    m %>% recover_types(d) %>%
      spread_draws(kappa[angle]) %>%
      mutate(angle = (angle-1)*90,
             angle = as_factor(angle)) %>%
      ungroup() %>%
      ggplot(aes(kappa, fill = angle)) + 
      geom_density(alpha = 0.5) + 
      geom_vline(data = tibble(angle = as_factor(0:1 * 90), kappa = gt$kappa), 
                 aes(xintercept = kappa, colour = angle), linetype = 1) + 
      ggthemes::scale_color_ptol() -> plt_abs_dir_kappa
    
    post_dir <- extract_post_fixed_abs_dir_comps(fit, d, n = 100) 
    
    map_df(unique(post_dir$.draw), compute_rho_phi) %>%
      ggplot(aes(x = phi, y = rho)) + geom_path(aes( group = .draw), alpha = 0.25) + 
      geom_path(data = ddir, aes(phi, z), colour = "darkred", linewidth = 2 ) -> plt_dir_dist
    
    plt <- plt + (plt_abs_dir_theta + plt_abs_dir_kappa + plt_dir_dist + plot_layout(widths = c(1,1,2)))
  }
  
  plt <- plt +
    plot_layout(guides = "collect")  & theme(legend.position = 'bottom')
  
  return(plt)
}

plot_model_indiv <- function(post, m, d) {
  
  
  
}
  
compute_rho_phi <- function(drw) {
  
  p <- filter(post_dir, .draw == drw)
  
  theta <- c(p$theta[1], p$theta[2], p$theta[1], p$theta[2])
  kappa <- c(p$kappa[1], p$kappa[2], p$kappa[1], p$kappa[2])
  phi = seq(1, 360)
  
  ddir <- tibble(phi = seq(1, 360),
                 .draw = drw,
                 rho = compute_all_von_mises(theta, kappa, phi))
  
  return(ddir)
}

plt_post_prior <- function(post, prior, var, xtitle, gt) {
  
  # function to plot the posterior against the prior. 
  # gt allows us to mark up the groundtruth (if available)
  
  prior_var = paste("prior", var, sep = "_")
  
  if (var == "p_floor") 
  {
    post %>% 
      ggplot() + 
      geom_rect(data = prior %>% 
                  median_hdci(exp(get(prior_var)), .width = c(0.53, 0.97)),
                aes(ymin = -Inf, ymax = Inf, xmin = .lower, xmax = .upper), 
                fill = "orange", alpha = 0.25) +  
      geom_density(aes(exp(get(var)), fill = condition), alpha = 0.5) +
      scale_x_continuous(xtitle) +
      coord_cartesian(xlim = c(0, 0.1)) -> plt
    
  } else {
    post %>% 
      ggplot() + 
      geom_rect(data = prior %>% 
                  median_hdci(get(prior_var), .width = c(0.53, 0.97)),
                aes(ymin = -Inf, ymax = Inf, xmin = .lower, xmax = .upper), 
                fill = "orange", alpha = 0.25) +  
      geom_density(aes(get(var), fill = condition), alpha = 0.5) +
      geom_vline(xintercept = gt, linetype = 1, colour = "red") +
      scale_x_continuous(xtitle) -> plt
    
  }
  
  return(plt)
  
}

plot_traceplots <- function(m) {
  
  bayesplot::mcmc_trace(m$draws(), pars = c("bA", "b_stick"))
  
}

plot_model_random <- function(post) {
  
  ggplot(post, aes(sig, fill = dist)) + 
    geom_density(alpha = 0.5) +
    facet_wrap(~param, scales = "free")
  
}

plt_beta_comp <- function(dat, var) {
  
  dat %>% filter(dim == var) %>%
    ggplot(aes(x, z, group = interaction(c, draw), colour = c)) +
    geom_path(alpha = 0.1) + 
    theme(legend.position = "none")  + 
    scale_y_continuous("likelihood")  -> plt
  
  if (var == "y") {
    plt <- plt + coord_flip() + 
      scale_x_reverse()  +
      xlab(var)
  }
  
  return(plt)
  
}

sample_beta <- function(.draw, dim, comp, a, b) {
  
  x <- seq(0.01, 0.99, 0.01)
  
  return(tibble(c = comp, 
                draw = .draw,
                dim = dim,
                x = x,
                z = dbeta(x, shape1 = a, shape2 = b)))
}

plot_init_sel <- function(post, d, pp=FALSE) {
  
  post$lambda %>%
    mutate(person = as_factor(person)) %>%
    ggplot(aes(lambda, person)) +
    stat_pointinterval(colour = "purple") -> plt_lambda
  
  if (pp==FALSE) {
  post$init_sel %>%
    sample_n(500) %>%
      pmap_df(sample_beta) %>%
      mutate(c = as_factor(c)) -> plt_dat
    
  } else {
    post$init_prior %>%
      sample_n(500) %>%
      mutate(a = a, 
             b = b) %>%
      pmap_df(sample_beta) %>%
      mutate(c = as_factor(c)) -> plt_dat
    
  }
   
    plt <- plt_lambda + plt_beta_comp(plt_dat, "y") + plt_beta_comp(plt_dat, "x")
   
    return(plt)
}
# 
# plot_model_spatial <- function(m, d) {
#   
#   m %>% 
#     spread_draws( sigma_dis[block], sigma_dir[block]) %>%
#     mutate(block = as_factor(block),
#            block = fct_recode(block, scarce = "1", equal = "2")) %>%
#     ungroup() -> post
#   
#   m %>% spread_draws(u[block, person],) %>%
#     mutate(param = block %% 4,
#            param = if_else(param == 0, 4, param),
#            block = (block / 4),
#            block = ceiling(block),
#            block = as_factor(block),
#            block = fct_recode(block, feature = "1", conjunction = "2"),
#            param = as_factor(param),
#            param = fct_recode(param, sigma_dir = "1", sigma_dis = "2")) %>%
#     rename(uz = "u") -> post_u
#   
#   full_join(post %>%
#               pivot_longer(c(sigma_dis, sigma_dir), names_to = "param", values_to = "u"), 
#             post_u, 
#             by = c("block", ".chain", ".iteration", ".draw", "param")) %>% 
#     mutate(uz = u + uz) %>%
#     group_by(person, block, param) %>%
#     summarise(uz = mean(uz), .groups = "drop") %>%
#     pivot_wider(names_from = "param", values_from = "uz") -> post_u
# 
#     
#   distances <- seq(0, 1, 0.025)
#   angles <- seq(0, 2*pi, 0.1)
#   
#   # draw distance tuning figure for foxed effect
#     pmap_dfr(select(post, block, sigma_dis), 
#              function(sigma_dis, t, block) tibble(t=t, block = block,
#                                                        q = exp(-sigma_dis * t)),
#              t = distances) %>%
#       group_by(block, t) %>%
#       median_hdci(q, .width = c(0.53, 0.97)) %>%
#       select(-q, -.point, -.interval) %>%
#       unite("interval", .lower, .upper) %>%
#       pivot_wider(names_from = ".width", values_from = "interval") %>%
#       separate(`0.53`, c("lower53", "upper53"), sep  = "_", convert = T)  %>%
#       separate(`0.97`, c("lower97", "upper97"), sep  = "_", convert = T) -> q
#     
#     # plot fixed effect distance -> weight function HDPI
#     q %>%
#       ggplot(aes(t)) +
#       geom_ribbon(aes(ymin = lower97, ymax = upper97, fill = block), alpha = 0.5) +
#       geom_ribbon(aes(ymin = lower53, ymax = upper53, fill = block), alpha = 0.7) +
#       scale_x_continuous("distance to target") +
#       scale_y_continuous("weighting") +
#       ggtitle("fixed effects") -> plt_dis
#     
#     pmap_dfr(select(post, block, direction_bias), 
#              function(direction_bias, a, block) tibble(a=a, block = block,
#                                                          q = (1 + direction_bias*cos(4*a))/(direction_bias+1)),
#              a = angles) %>%
#       group_by(block, a) %>%
#       median_hdci(q, .width = c(0.53, 0.97)) %>%
#       select(-q, -.point, -.interval) %>%
#       unite("interval", .lower, .upper) %>%
#       pivot_wider(names_from = ".width", values_from = "interval") %>%
#       separate(`0.53`, c("lower53", "upper53"), sep  = "_", convert = T)  %>%
#       separate(`0.97`, c("lower97", "upper97"), sep  = "_", convert = T) -> q_direction
#     
#     # plot fixed effect distance -> weight function HDPI
#     q_direction %>%
#       ggplot(aes(a)) +
#       geom_ribbon(aes(ymin = lower97, ymax = upper97, fill = block), alpha = 0.5) +
#       geom_ribbon(aes(ymin = lower53, ymax = upper53, fill = block), alpha = 0.7) +
#       scale_x_continuous("distance to target") +
#       scale_y_continuous("weighting") +
#       ggtitle("fixed effects") -> plt_dir
#     
#     # now plot indiv prox fall offs:
#     pmap_dfr(select(post_u, block, sigma_dis, person), 
#              function(sigma_dis, t, block, person) tibble(t=t, block = block, person = person,
#                                                          q = exp(-sigma_dis * t)),
#              t = distances) %>%
#       ggplot(aes(x = t, y = q, colour = block, group = interaction(block, person))) + 
#       geom_path(alpha = 0.5) +
#       scale_x_continuous("distance to target") +
#       scale_y_continuous("weighting") +
#       ggtitle("indiv. differences") -> plt_dis_p
#     
#     
#     
#      plt_dis + plt_dis_p + plt_dir
#   
# }


