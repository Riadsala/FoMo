# compile all of our models in one place.
# this is useful for checking everything works!

library(cmdstanr)


cmdstan_model("simple/fomo1_0.stan", force = T)
cmdstan_model("multi_level/fomo1_0.stan", force = T)
cmdstan_model("simulate/fomo1_0.stan", force = T)

cmdstan_model("simple/fomo1_2.stan", force = T)
cmdstan_model("multi_level/fomo1_2.stan", force = T)
cmdstan_model("simulate/fomo1_2.stan", force = T)

cmdstan_model("simple/fomo1_3.stan", force = T)
cmdstan_model("multi_level/fomo1_3.stan", force = T)
cmdstan_model("simulate/fomo1_3.stan", force = T)

cmdstan_model("multi_level/fomo1_4.stan", force = T)
cmdstan_model("simulate/fomo1_4.stan", force = T)