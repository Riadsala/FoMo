/* FoMo V1.3 - multi-level

Adds in rho_phi - abs direction


*/

functions {

  #include /../include/FoMo_functions.stan

  vector compute_weights(
    real u_a, real u_s, real u_delta, real u_psi, vector log_theta, real kappa,
    vector item_class, vector match_prev_item, vector delta, vector psi, vector phi, 
    int n, int n_targets, vector remaining_items, real os) {

    vector[n_targets] weights;
    
    // set the weight of each target to be its class weight
    weights = log_inv_logit(u_a * to_vector(item_class));

    // multiply weights by stick/switch preference
    weights += log_inv_logit(u_s * match_prev_item); 

    // calculate by spatial weights
    weights += compute_spatial_weights(n, n_targets, 
      u_delta, u_psi, log_theta, kappa,
      delta, psi, phi, os);
        
    // remove already-selected items, and standarise to sum = 1 
    weights = standarise_weights(exp(weights), n_targets, remaining_items); 

    return(weights);

  }

  vector compute_spatial_weights(int n, int n_targets, 
    real rho_delta,real rho_psi, vector log_theta, real kappa,
    vector delta, vector psi, vector phi, real os) {

    // computes spatial weights
    // for FoMo1.3, this includes proximity and absolute direction
    vector[n_targets] prox_weights, reldir_weights, absdir_weights;
  
    // apply spatial weighting
    prox_weights   = compute_prox_weights(n, n_targets, 
                                 rho_delta, delta);

    reldir_weights = compute_reldir_weights(n, n_targets, 
                                 rho_psi, psi);

    absdir_weights = compute_absdir_weights_fixed_kappa(n, n_targets, 
                                 log_theta, kappa, phi);

    // return the dot product of the weights
    return(prox_weights + absdir_weights + reldir_weights);

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

  // theta is a 4D vector containing the mixture weights for our direction model
  array[K, 4] real log_theta; // mixing values for abs directions

  ///////////////////////////////
  // random effects
  ///////////////////////////////

  // random effect variances
  vector<lower=0>[4*K] sigma_u;
  vector<lower=0>[4*K] sigma_w;
  // declare L_u to be the Choleski factor of a 4*K x 4*K correlation matrix
  cholesky_factor_corr[4*K] L_u;
  //cholesky_factor_corr[4*K] L_w; // skip var-cov matrix for now

  // random effect matrix
  matrix[4*K, L] z_u; // 4*K as we have three params
  matrix[4*K, L] z_w; // 4*K as we have four directions
  
}

transformed parameters {

  /* 
  combine fixed and random effects
  we do this here so that the code in the model{} block
  is easier to read
  */

  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[4*K, L] u = diag_pre_multiply(sigma_u, L_u) * z_u;
  //matrix[4*K, L] w = diag_pre_multiply(sigma_w, L_w) * z_w;

  // add fixed and random effects together
  // create empty arrays for everything
  array[K] vector[L] u_a, u_stick, u_psi, u_delta;

  for (k in 1:K) {
    u_a[k]     = to_vector(b_a[k]       + u[4*(k-1)+1]);
    u_stick[k] = to_vector(b_stick[k]   + u[4*(k-1)+2]);
    u_delta[k] = to_vector(rho_delta[k] + u[4*(k-1)+3]);
    u_psi[k]   = to_vector(rho_psi[k] + u[4*(k-1)+4]);
  }

  // now work out thet u_log_theta
  array[K, L] vector[4] u_log_theta;

  for (k in 1:K) {
    for (l in 1:L) {

      u_log_theta[k, l, 1] = (log_theta[k, 1] + z_w[4*(k-1)+1, l] .* sigma_w[4*(k-1)+1]);
      u_log_theta[k, l, 2] = (log_theta[k, 2] + z_w[4*(k-1)+2, l] .* sigma_w[4*(k-1)+2]);
      u_log_theta[k, l, 3] = (log_theta[k, 3] + z_w[4*(k-1)+3, l] .* sigma_w[4*(k-1)+3]);
      u_log_theta[k, l, 4] = (log_theta[k, 4] + z_w[4*(k-1)+4, l] .* sigma_w[4*(k-1)+4]);

    }  
  }
}

model {

  /////////////////////////////////////////////////////
  // Define Priors
  ////////////////////////////////////////////////////

  // priors for random effects
  sigma_u ~ exponential(1);
  sigma_w ~ exponential(1);

  // LKJ prior for the correlation matrix
  L_u~ lkj_corr_cholesky(1.5); 
  //L_w ~ lkj_corr_cholesky(1.5);

  to_vector(z_u) ~ normal(0,1);
  to_vector(z_w) ~ normal(0,1);

  for (ii in 1:K) {
    // priors for fixed effects
    target += normal_lpdf(b_a[ii]       | 0, prior_sd_b_a);
    target += normal_lpdf(b_stick[ii]   | 0, prior_sd_b_stick);
    target += normal_lpdf(rho_delta[ii] | prior_mu_rho_delta, prior_sd_rho_delta);
    target += normal_lpdf(rho_psi[ii]   | prior_mu_rho_psi, prior_sd_rho_psi);
    target += normal_lpdf(log_theta[ii] | 0, 2);
  }

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  vector[n_targets] weights;

  // some IDs for trial, participant, and condition
  int t, z, x; 

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  for (ii in 1:N) {

    t = trial[ii];
    z = Z[t];
    x = X[t];
 
    weights = compute_weights(
      u_a[x, z], u_stick[x, z], u_delta[x, z], u_psi[x, z], u_log_theta[x, z], kappa,
      to_vector(item_class[t]), S[ii], delta[ii], psi[ii], phi[ii],
      found_order[ii], n_targets, remaining_items[ii],
      grid_offset[x]); 

    target += log(weights[Y[ii]]);
  
  }
}

generated quantities {
  // here we  can output our prior distritions
  real prior_b_a = normal_rng(prior_mu_b_a, prior_sd_b_a);
  real prior_b_stick = normal_rng(prior_mu_b_stick, prior_sd_b_stick);
  real prior_rho_delta = normal_rng(prior_mu_rho_delta, prior_sd_rho_delta);
  real prior_rho_psi = normal_rng(prior_mu_rho_psi, prior_sd_rho_psi);

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
        u_a[x, z], u_stick[x, z], u_delta[x, z], u_psi[x, z], u_log_theta[x, z], kappa,
        to_vector(item_class[t]), S[ii], delta[ii], psi[ii], phi[ii],
        found_order[ii], n_targets, remaining_items[ii],
      grid_offset[x]); 

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
    vector[n_targets] S_j, psi_j, phi_j, delta_j;
    vector[n_targets] weights;  // class weight for each target

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
        phi_j = rep_vector(1, n_targets);
        psi_j = rep_vector(1, n_targets);

            
        if (ii > 1) {

          S_j     = compute_matching(item_class[t], n_targets, Q[t, ], ii);
          delta_j = compute_prox(item_x[t], item_y[t], n_targets, Q[t, ], ii);
          psi_j   = compute_reldir(item_x[t], item_y[t], n_targets, Q[t, ], ii);
          phi_j   = compute_absdir(item_x[t], item_y[t], n_targets, Q[t, ], ii); 
        }

        weights = compute_weights(
          u_a[x, z], u_stick[x, z], u_delta[x, z], u_psi[x, z], u_log_theta[x, z], kappa,
          to_vector(item_class[t]), S_j, delta_j, psi_j, phi_j,
          found_order[ii], n_targets, remaining_items_j,
          grid_offset[x]); 

        Q[t, ii] = categorical_rng(weights);

        // update remaining_items2
        remaining_items_j[Q[t, ii]] = 0;
 
      }   
    }
  } 
}