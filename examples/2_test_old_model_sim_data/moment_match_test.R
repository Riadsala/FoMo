library(tidyverse)
library(cmdstanr)

mod <- cmdstan_model("test_stan_code.stan", force_recompile = TRUE, cpp_options = list(stan_threads = TRUE))

iter = 100

d = list(y = c(4,5,6,7,8),
         N = 5)

m <- mod$sample(data = d, 
                           chains = 4, parallel_chains = 4, threads = 4,
                           refresh = 0, 
                           iter_warmup = iter, iter_sampling = iter,
                           sig_figs = 3)


loo_test <- m$loo(moment_match = TRUE)
