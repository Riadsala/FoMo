/* FoMo V1.2 - multi-level

Removes relative direction:

b_a, b_stick, rho_delta

*/

functions {

  #include /../include/FoMo_functions.stan

  vector compute_weights(
    real u_a, real u_s, real u_delta,
    vector item_class, vector match_prev_item, vector delta,
    int n, int n_targets, vector remaining_items) {

    vector[n_targets] weights;
    
    // set the weight of each target to be its class weight
    weights = log_inv_logit(u_a * to_vector(item_class));

    // multiply weights by stick/switch preference
    weights += log_inv_logit(u_s * match_prev_item); 

    // calculate by spatial weights
    weights += compute_spatial_weights(
      n, n_targets, u_delta, delta);

    // remove already-selected items, and standarise to sum = 1 
    weights = standarise_weights(exp(weights), n_targets, remaining_items); 

    return(weights);

  }

  vector compute_spatial_weights(int n, int n_targets, 
    real rho_delta, vector delta) {

    // computes spatial weights
    // for FoMo1.0, this includes proximity and relative direction
    vector[n_targets] prox_weights;


    // apply spatial weighting
    prox_weights   = compute_prox_weights(n, n_targets, 
     rho_delta, delta);

    // return the dot product of the weights
    return(prox_weights);

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

  ///////////////////////////////
  // random effects
  ///////////////////////////////
  // random effect variances
  vector<lower=0>[3*K] sigma_u;
  // declare L_u to be the Choleski factor of a 3*Kx3*K correlation matrix
  cholesky_factor_corr[3*K] L_u;
  // random effect matrix
  matrix[3*K,L] z_u; 
  
}

transformed parameters {

  /* 
  combine fixed and random effects
  we do this here so that the code in the model{} block
  is easier to read
  */

  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[3*K, L] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u;

  // add fixed and random effects together
  // create empty arrays for everything
  array[K] vector[L] u_a, u_stick, u_delta, u_psi;
  // create!
  for (kk in 1:K) {
    u_a[kk]     = to_vector(b_a[kk]       + u[3*(kk-1)+1]);
    u_stick[kk] = to_vector(b_stick[kk]   + u[3*(kk-1)+2]);
    u_delta[kk] = to_vector(rho_delta[kk] + u[3*(kk-1)+3]);
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
  for (ii in 1:N) {

    t = trial[ii];
    z = Z[t];
    x = X[t];

    weights = compute_weights(
      u_a[x, z], u_stick[x, z], u_delta[x, z],
      to_vector(item_class[t]), S[ii], delta[ii],
      found_order[ii], n_targets, remaining_items[ii]); 

    target += log(weights[Y[ii]]);

  }
}

generated quantities {
  // here we  can output our prior distritions
  real prior_b_a = normal_rng(prior_mu_b_a, prior_sd_b_a);
  real prior_b_stick = normal_rng(prior_mu_b_stick, prior_sd_b_stick);
  real prior_rho_delta = normal_rng(prior_mu_rho_delta, prior_sd_rho_delta);

  array[N] int P;
  array[N] real log_lik;

  // for trial level predictions, we have to remember that we do not have a stopping rule yet
  // so we will simply collect all of the targets
  array[n_trials, n_targets] int Q;

  //////////////////////////////////////////////////////////////////////////////
  /* 
  i) Step through data item-selection at a time and compute which item the 
  model selects next. This allows us to compute accuracy
  */
  //////////////////////////////////////////////////////////////////////////////
  {
    // some counters and index variables, etc.
    vector[n_targets] weights;  // class weight for teach target
    
    // some IDs for trial, condition, and condition
    int t, z, x; 

    //////////////////////////////////////////////////
    // // step through data row by row and define LLH
    //////////////////////////////////////////////////
    for (ii in 1:N) {

      t = trial[ii];
      z = Z[t];
      x = X[t];

      weights = compute_weights(
        u_a[x, z], u_stick[x, z], u_delta[x, z], 
        to_vector(item_class[t]), S[ii], delta[ii],
        found_order[ii], n_targets, remaining_items[ii]); 

      P[ii] = categorical_rng(weights);
      log_lik[ii] = log(weights[Y[ii]]);

    }
  }

  //////////////////////////////////////////////////////////////////////////////
  /* 
  ii) Simulate each trial, start to finish. This allows us to compuare the model
  to human participants in terms of run statistics and selection vectors
  */
  //////////////////////////////////////////////////////////////////////////////
  {
    vector[n_targets] remaining_items_j;
    vector[n_targets] S_j, psi_j, delta_j;
    vector[n_targets] weights;  // 

    //for each trial
    for (t in 1:n_trials) {

      int z = Z[t];
      int x = X[t];

      // first, set things up!
      remaining_items_j = rep_vector(1, n_targets);


      // simulate a trial!
      for (ii in 1:n_targets) {

        // delta_j: distance to previously selected item
        S_j = rep_vector(0, n_targets);
        delta_j = rep_vector(1, n_targets);

        if (ii > 1) 
        {
          S_j     = compute_matching(item_class[t], n_targets, Q[t, ], ii);
          delta_j = compute_prox(item_x[t], item_y[t], n_targets, Q[t, ], ii);
        }

        weights = compute_weights(
          u_a[x, z], u_stick[x, z], u_delta[x, z],
          to_vector(item_class[t]), S_j, delta_j,
          found_order[ii], n_targets, remaining_items_j); 

        Q[t, ii] = categorical_rng(weights);

        // update remaining_items2
        remaining_items_j[Q[t, ii]] = 0;

      }
    }

  }
} 