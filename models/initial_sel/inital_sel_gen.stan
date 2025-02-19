/* FoMo V1.0 - multi-level

Includes the core parameters:

b_a, b_stick, rho_delta, rho_psi

*/

functions {

  vector beta_weight(vector X, real a, real b) {

    int n = size(X);
    vector[n] Z;
    Z = X^(a-1).*(1-X^(b-1));
    //Z = Z / beta(a, b); 
    return(Z);
  }

  #include /../include/FoMo_functions.stan

  vector compute_weights(
    real u_a, real u_s, real u_delta, real u_psi,
    vector x, vector y,
    real a_x, real a_y, real b_x, real b_y, real c_x, real c_y, real lambda,
    vector item_class, vector match_prev_item, vector delta, vector psi,
    int n, int n_targets, vector remaining_items) {

    vector[n_targets] weights;
    
    // set the weight of each target to be its class weight
    weights = log_inv_logit(u_a * to_vector(item_class));

    // multiply weights by stick/switch preference
    weights += log_inv_logit(u_s * match_prev_item); 

    if (n == 1) {

      vector[n_targets] w1, w2;

      // compute intial selection weights
      w1 = beta_weight(x, (c_x), (c_x)) .* beta_weight(y, (c_y), (c_y));
      w2 = beta_weight(x, (a_x), (b_x)) .* beta_weight(y, (a_y), (b_y));

      weights += inv_logit(lambda) * w1 + (1-inv_logit(lambda)) * w2;
      /*
      print("lambda: ", lambda, " ----- logit lambda: ", logit(lambda));
      print("w1 and w2");
      print(w1);
      print(w2);
      print(weights);
      */

      // remove already-selected items, and standarise to sum = 1 
      weights = standarise_weights(weights, n_targets, remaining_items); 



    } else {

      // calculate by spatial weights
      weights += compute_spatial_weights(
        n, n_targets, u_delta, u_psi, delta, psi);

      // remove already-selected items, and standarise to sum = 1 
      weights = standarise_weights(exp(weights), n_targets, remaining_items); 
  }
        
    

    return(weights);

  }

  vector compute_spatial_weights(int n, int n_targets, 
    real rho_delta, real rho_psi, vector delta, vector psi) {

    // computes spatial weights
    // for FoMo1.0, this includes proximity and relative direction
    vector[n_targets] prox_weights;
    vector[n_targets] reldir_weights;

    // apply spatial weighting
    prox_weights   = compute_prox_weights(n, n_targets, 
                                 rho_delta, delta);
    
    reldir_weights = compute_reldir_weights(n, n_targets, 
                                 rho_psi, psi);

    // return the dot product of the weights
    return(prox_weights + reldir_weights);

  }
}


data {
  int <lower = 1> N; // total number of selected targets over the whole experiment
  int <lower = 1> L; // number of participant levels 
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
  array[n_trials] int <lower = 1, upper = L> Z; // random effect levels

  array[n_trials, n_targets] int <lower = -1, upper = 1> item_class; // target class, one row per trial
  array[N] vector<lower = -1, upper = 1>[n_targets] S; // stick/switch (does this targ match prev targ) 
  
  array[N] int<lower = 1, upper = n_trials> trial; // what trial are we on? 

  // read in priors
  real prior_mu_b_a; // param for class weight prior
  real prior_sd_b_a; // param for class weight prior
  real prior_mu_b_stick; // prior for sd for bS
  real prior_sd_b_stick; // prior for sd for bS
  real prior_mu_rho_delta;
  real prior_sd_rho_delta;
  real prior_mu_rho_psi;
  real prior_sd_rho_psi;

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

  /* in order to allow for correlations between the
  variables, these are all stored in a list
  these include b_a, bS (stick weight), and the two spatial 
  sigmas, along with the floor (chance of selectin an 
  item at random)
  */
  array[K] real b_a; // weights for class A compared to B  
  array[K] real b_stick; // stick-switch rates 
  array[K] real<lower = 0> rho_delta; // distance tuning
  array[K] real rho_psi; // direction tuning

  ///////////////////////////////
  // random effects
  ///////////////////////////////

  // random effect variances
  vector<lower=0>[4*K] sigma_u;
  // declare L_u to be the Choleski factor of a 3*Kx3*K correlation matrix
  cholesky_factor_corr[4*K] L_u;
  // random effect matrix
  matrix[4*K,L] z_u; 

  ///////////////////////////////
  // initial item selection
  ///////////////////////////////
  array[K, L] real lambda; // mixture weight for initial strategy 
  real<lower=0> a_x;
  real<lower=0> a_y;
  real<lower=0> b_x;
  real<lower=0> b_y;
  real<lower=0> c_x;
  real<lower=0> c_y;

}

transformed parameters {

  /* 
  combine fixed and random effects
  we do this here so that the code in the model{} block
  is easier to read
  */

  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[4*K, L] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u;

  // add fixed and random effects together
  // create empty arrays for everything
  array[K] vector[L] u_a, u_stick, u_delta, u_psi;
  // create!
  for (kk in 1:K) {
    u_a[kk]     = to_vector(b_a[kk]       + u[4*(kk-1)+1]);
    u_stick[kk] = to_vector(b_stick[kk]   + u[4*(kk-1)+2]);
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
  L_u ~ lkj_corr_cholesky(1.5); // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0,1);

  for (ii in 1:K) {
    // priors for fixed effects
    target += normal_lpdf(b_a[ii]       | 0, prior_sd_b_a);
    target += normal_lpdf(b_stick[ii]   | 0, prior_sd_b_stick);
    target += normal_lpdf(rho_delta[ii] | prior_mu_rho_delta, prior_sd_rho_delta);
    target += normal_lpdf(rho_psi[ii]   | prior_mu_rho_psi, prior_sd_rho_psi);
  }


  // priors for inital item selection
  a_x ~ normal(0, 1);
  a_y ~ normal(0, 1);

  b_x ~ normal(0, 1);
  b_y ~ normal(0, 1);

  c_x ~ normal(0, 1);
  c_y ~ normal(0, 1);

  for (k in 1:2) {
    lambda[k] ~ normal(0, 1);
  }

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  vector[n_targets] weights;

  // some IDs for trial, condition, and condition
  int t, z, x; 

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  ////////////////////////////////////////////////// 
  //print(a_x, "-", b_x, "-", b_x, "-", b_y, "-", c_x, "-", c_y); 
  for (ii in 1:N) {

    
    
    //print(ii);

    t = trial[ii];
    z = Z[t];
    x = X[t];

    //print(ii, " - ", item_x[t]);
 
    weights = compute_weights(
      u_a[x, z], u_stick[x, z], u_delta[x, z], u_psi[x, z],
      item_x[t], item_y[t],
      a_x, a_y, b_x, b_y, c_x, c_y, lambda[x, z],
      to_vector(item_class[t]), S[ii], delta[ii], psi[ii],
      found_order[ii], n_targets, remaining_items[ii]); 

    //print(weights);

    target += log(weights[Y[ii]]);
   
  }
}

generated quantities {


}