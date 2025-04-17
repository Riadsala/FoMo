/* 

FoMo V1.0 (multi-level)

This is the new implementation of the model from 
Clarke et al (2022), Comp Bio.

Includes the core parameters:
b_a, b_s, rho_delta, rho_psi

*/

functions {

  #include /../include/FoMo_functions.stan

}

data {
  int <lower = 1> N; // total number of selected targets over the whole experiment
  int <lower = 1> K; // number of experimental conditions  
  int <lower = 1> L; // number of participant levels 
  
  int <lower = 1> n_trials;  // total number of trials (overall)
  int <lower = 1> n_classes; // number of target classes - we assume this is constant over n_trials
  int <lower = 1> n_targets; // total number of targets per trial

  array[N] int <lower = 1> Y; // target IDs - which target was selected here? This is what we predict
  array[n_trials] int <lower = 1, upper = K> X; // trial features (ie, which condition are we in)
  array[n_trials] int <lower = 1, upper = L> Z; // random effect levels

  array[N] int<lower = 1, upper = n_trials> trial; // what trial are we on? 
  array[N] int <lower = 0, upper = n_targets> found_order; // = 1 if starting a new trial, 0 otherwise
  
  // (x, y) coordinates of each target
  array[n_trials] vector<lower=0,upper=1>[n_targets] item_x;
  array[n_trials] vector<lower=0,upper=1>[n_targets] item_y;

  // item class features
  array[n_trials, n_targets] int <lower = -1, upper = 1> item_class; // target class, one row per trial
  array[N] vector<lower = -1, upper = 1>[n_targets] S; // stick/switch (does this targ match prev targ) 

  // pre-computed inter-item features
  array[N] vector<lower = 0>[n_targets] delta; // distance measures
  array[N] vector[n_targets] psi; // direction measures (relative)
  array[N] vector[n_targets] phi; // direction measures (absolute)

  // read in priors
  // suggested values given in comments
  real prior_mu_b_a; // = 0, prior for salience of item class A compared to B
  real prior_sd_b_a; // = 1, uncertainty for b_a prior
  real prior_mu_b_s; // = 0, prior for b_s, item class sticking v switching
  real prior_sd_b_s; // = 1, uncertainty for b_s prior
  real prior_mu_rho_delta; // = 15, negexp fall off due to proximity
  real prior_sd_rho_delta; // = 5, uncertainty around rho_delta
  real prior_mu_rho_psi; // = 0, "momentum"
  real prior_sd_rho_psi; // = 0.5, uncertainty around rho_psi

}

transformed data {

  // pre-compute array giving remaining items for each selection
  array[N] vector[n_targets] remaining_items;
  remaining_items = calc_remaining_items(N, n_targets, Y, found_order);

}

parameters {

  // these are the parameters we wish to fit to the data

  ////////////////////////////////////
  // fixed effects
  ////////////////////////////////////
  array[K] real b_a; // weights for class A compared to B  
  array[K] real b_s; // stick-switch rates 
  array[K] real<lower = 0> rho_delta; // distance tuning
  array[K] real rho_psi; // direction tuning

  ///////////////////////////////
  // random effects
  ///////////////////////////////
  // random effect variances: 
  // 4*K as we have four fixed effect parameters x K conditions
  vector<lower=0>[4*K] sigma_u;
  cholesky_factor_corr[4*K] L_u;
  // random effect matrix
  matrix[4*K,L] z_u; 
  
}

transformed parameters {

  /* 
  combine fixed and random effects
  we do this here so that the code in the model{} block is easier to read
  */

  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[4*K, L] u = diag_pre_multiply(sigma_u, L_u) * z_u;

  // create empty arrays for everything
  array[K] vector[L] u_a, u_stick, u_delta, u_psi;
  // calculate
  for (kk in 1:K) {
    u_a[kk]     = to_vector(b_a[kk]       + u[4*(kk-1)+1]);
    u_stick[kk] = to_vector(b_s[kk]       + u[4*(kk-1)+2]);
    u_delta[kk] = to_vector(rho_delta[kk] + u[4*(kk-1)+3]);
    u_psi[kk]   = to_vector(rho_psi[kk]   + u[4*(kk-1)+4]);
  }
}

model {

  /////////////////////////////////////////////////////
  // Define Priors
  ////////////////////////////////////////////////////

  // priors for random effects
  sigma_u ~ exponential(1);
  L_u ~ lkj_corr_cholesky(2); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0, 1); // centred prior for random effects, so this should always be N(0,1)

  // priors for fixed effects
  for (kk in 1:K) {
    target += normal_lpdf(b_a[kk]       | prior_mu_b_a, prior_sd_b_a);
    target += normal_lpdf(b_s[kk]       | prior_mu_b_s, prior_sd_b_s);
    target += normal_lpdf(rho_delta[kk] | prior_mu_rho_delta, prior_sd_rho_delta);
    target += normal_lpdf(rho_psi[kk]   | prior_mu_rho_psi, prior_sd_rho_psi);
  }

  // create some variables
  vector[n_targets] weights;
  int t, z, x; // trial, person, and condition

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  for (ii in 1:N) {

    t = trial[ii];
    z = Z[t];
    x = X[t];
 
    weights = compute_weights_v10(
      u_a[x, z], u_stick[x, z], u_delta[x, z], u_psi[x, z],
      to_vector(item_class[t]), S[ii], delta[ii], psi[ii],
      found_order[ii], n_targets, remaining_items[ii]); 

    // get likelihood of item selection
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
