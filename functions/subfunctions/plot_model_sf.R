# the functions here are not designed to be called directly.
# they're all used by the tools in ../plot_model.R

##################################################################################
# subfunctions for plot_model_fixed()

plot_cts_params <- function(post,  gt=NULL, clist=NULL) {
  
  # this function is currently unused?????
  
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

plot_model_weights <- function(post, params) {
  
  plts <- map(params, create_weight_plot, post = post)
  
  plt <- wrap_plots(plts) + plot_layout(guides = "collect")
  
  return(plt)
  
}


