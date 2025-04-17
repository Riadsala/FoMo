/* 

FoMo V1.3 (single-level)

This model adds absolute direction (psi)

Includes the core parameters:
b_a, b_s, rho_delta, rho_psi and 
a set of theta mixture weights

kappa is passed in as a hyper parameter

*/


functions {

  #include /../include/FoMo_functions.stan

}

data {
  int <lower = 1> N; // total number of selected targets over the whole experiment
  int <lower = 1> K; // number of experimental conditions  

  int <lower = 1> n_trials;  // total number of trials (overall)
  int <lower = 1> n_classes; // number of target classes - we assume this is constant over n_trials
  int <lower = 1> n_targets; // total number of targets per trial
  array[N] int <lower = 0, upper = n_targets> found_order; // = 1 is starting a new trial, 0 otherwise

  array[N] int <lower = 1> Y; // target IDs - which target was selected here? This is what we predict

  // (x, y) coordinates of each target
  array[n_trials] vector<lower=0,upper=1>[n_targets] item_x;
  array[n_trials] vector<lower=0,upper=1>[n_targets] item_y;

  array[N] vector<lower = 0>[n_targets] delta; // distance measures
  array[N] vector[n_targets] psi; // direction measures (relative)
  array[N] vector[n_targets] phi; // direction measures (absolute)

  array[n_trials] int <lower = 1, upper = K> X; // trial features (ie, which condition are we in)

  array[n_trials, n_targets] int <lower = -1, upper = 1> item_class; // target class, one row per trial
  array[N] vector<lower = -1, upper = 1>[n_targets] S; // stick/switch (does this targ match prev targ) 
  
  array[N] int<lower = 1, upper = n_trials> trial; // what trial are we on? 

  // read in priors
  real prior_mu_b_a; // param for class weight prior
  real prior_sd_b_a; // param for class weight prior
  real prior_mu_b_s; // prior for sd for bS
  real prior_sd_b_s; // prior for sd for bS
  real prior_mu_rho_delta;
  real prior_sd_rho_delta;
  real prior_mu_rho_psi;
  real prior_sd_rho_psi;
  real prior_theta_lambda;

  // pass in kappa hyper-parameter
  real<lower = 0> kappa;
  array[K] real grid_offset; // angular grid offset (ie, should we rotate our coordiante acis)
  
}

transformed data{

  array[N] vector[n_targets] remaining_items;
  remaining_items = calc_remaining_items(N, n_targets, Y, found_order);

}

parameters {
  // These are all the parameters we want to fit to the data

  ////////////////////////////////////
  // fixed effects
  ////////////////////////////////////

  array[K] real b_a; // weights for class A compared to B  
  array[K] real b_s; // stick-switch rates 
  array[K] real<lower = 0> rho_delta; // distance tuning
  array[K] real rho_psi; // direction tuning

   // theta is a 4D vector containing the mixture weights for our direction model
  array[K] vector[4] log_theta; // mixing proportions for abs directions
}

model {

  /////////////////////////////////////////////////////
  // Define Priors
  ////////////////////////////////////////////////////

  for (kk in 1:K) {
    // priors for fixed effects
    target += normal_lpdf(b_a[kk]       | 0, prior_sd_b_a);
    target += normal_lpdf(b_s[kk]   | 0, prior_sd_b_s);
    target += normal_lpdf(rho_delta[kk] | prior_mu_rho_delta, prior_sd_rho_delta);
    target += normal_lpdf(rho_psi[kk]   | prior_mu_rho_psi, prior_sd_rho_psi);
    target += normal_lpdf(log_theta[kk] | 0, 2);
  }

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  vector[n_targets] weights;

  // some IDs for trial, condition, and condition
  int t, x; 

  for (ii in 1:N) {

    t = trial[ii];
    x = X[t];
 
    weights = compute_weights_v13(
      b_a[x], b_s[x], rho_delta[x], rho_psi[x], log_theta[x], kappa,
      to_vector(item_class[t]), S[ii], delta[ii], psi[ii], phi[ii],
      found_order[ii], n_targets, remaining_items[ii], grid_offset[1]); 

    target += weights[Y[ii]];
   
  }
}

generated quantities {
  // here we  can output our prior distritions
  real prior_b_a = normal_rng(prior_mu_b_a, prior_sd_b_a);
  real prior_b_s = normal_rng(prior_mu_b_s, prior_sd_b_s);
  real prior_rho_delta = normal_rng(prior_mu_rho_delta, prior_sd_rho_delta);
  real prior_rho_psi = normal_rng(prior_mu_rho_psi, prior_sd_rho_psi);

}