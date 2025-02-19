/* FoMo V1.0 - single-level

Includes the core parameters:

b_a, b_stick, rho_delta, rho_psi

*/

functions {

  #include /../include/FoMo_functions.stan

  vector compute_weights(
    real b_a, real b_s, real rho_delta, real rho_psi,
    vector item_class, vector match_prev_item, vector delta, vector psi,
    int n, int n_targets, vector remaining_items) {

    vector[n_targets] weights;
    
    // set the weight of each target to be its class weight
    weights = log_inv_logit(b_a * to_vector(item_class));

    // multiply weights by stick/switch preference
    weights += log_inv_logit(b_s * match_prev_item); 

    // calculate by spatial weights
    weights += compute_spatial_weights(
      n, n_targets, rho_delta, rho_psi, delta, psi);
        
    // remove already-selected items, and standarise to sum = 1 
    weights = standarise_weights(exp(weights), n_targets, remaining_items); 

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

  int <lower = 1> n_targets; // total number of targets per trial
  // (x, y) coordinates of each target
  vector<lower=0,upper=1>[n_targets] item_x;
  vector<lower=0,upper=1>[n_targets] item_y;

  array[N] int <lower = 1, upper = K> X; // trial features (ie, which condition are we in)
  array[N] int <lower = 1, upper = L> Z; // random effect levels
}

transformed data{



}

parameters {
  // These are all the parameters we want to fit to the data

  ////////////////////////////////////
  // fixed effects
  ////////////////////////////////////

  array[K, L] real lambda; // mixture weight for initial strategy 
  real<lower=0> a_x;
  real<lower=0> a_y;
  real<lower=0> b_x;
  real<lower=0> b_y;
  real<lower=0> c_x;
  real<lower=0> c_y;
  
}

model {

  // priors for mean (x, y)
  a_x ~ normal(1, 1);
  a_y ~ normal(1, 1);

  b_x ~ normal(5, 2);
  b_y ~ normal(5, 2);

  c_x ~ normal(2, 1);
  c_y ~ normal(2, 1);

  for (k in 1:2) {
    lambda[k] ~ normal(0, 1);
  }

  for (n in 1:N) {
    
  int k = X[n];
  int z = Z[n];
    
  target += log_mix(inv_logit(lambda[k, z]),
    beta_lpdf(item_x[n] | c_x, c_x),
    beta_lpdf(item_x[n] | a_x, b_x)); 

  target += log_mix(inv_logit(lambda[k, z]),
    beta_lpdf(item_y[n] | c_y, c_y),
    beta_lpdf(item_y[n] | a_y, b_y)); 

  }

 }
