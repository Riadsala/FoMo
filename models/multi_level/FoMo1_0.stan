// spatial foraging project

functions {

  #include /../include/FoMo_functions.stan

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
    return(prox_weights .* reldir_weights);

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
  matrix<lower = -1, upper = 1>[n_trials, n_targets] item_class; // target class, one row per trial
  array[N] vector<lower = -1, upper = 1>[n_targets] S; // stick/switch (does this targ match prev targ) 
  array[N] int <lower = 1, upper = L> Z; // random effect levels
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

  // parameters for simulation (generated quantities)
  int<lower = 0> n_trials_to_sim;
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

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  vector[n_targets] weights;
  vector[n_targets] spatial_weights;

  // some counters and index variables, etc.
  int t; // trial counter
  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  for (ii in 1:N) {

    t = trial[ii];
 
    // set the weight of each target to be its class weight
    weights = (u_a[X[t], Z[ii]]) * to_vector(item_class[t]) ;

    // multiply weights by stick/switch preference
    weights = inv_logit(weights) .* inv_logit(u_stick[X[t], Z[ii]] * S[ii]); 

    spatial_weights = compute_spatial_weights(found_order[ii], n_targets,
       u_delta[X[t], Z[ii]], u_psi[X[t], Z[ii]], delta[ii], psi[ii]);

    weights = weights .* spatial_weights;
        
    // remove already-selected items, and standarise to sum = 1 
    weights = standarise_weights(weights, n_targets, remaining_items[ii]);   
    //print("Y ", Y[ii]);
    //print("item y ", item_y[t]);
    //print(weights);

    target += log((weights)[Y[ii]]);

   
  }
}

generated quantities {
  // here we  can output our prior distritions
  real prior_b_a = normal_rng(prior_mu_b_a, prior_sd_b_a);
  real prior_b_stick = normal_rng(prior_mu_b_stick, prior_sd_b_stick);
  real prior_rho_delta = normal_rng(prior_mu_rho_delta, prior_sd_rho_delta);
  real prior_rho_psi = normal_rng(prior_mu_rho_psi, prior_sd_rho_psi);


  array[N] int P;
  array[N] real W;

  // for trial level predictions, we have to remember that we do not have a stopping rule yet
  // so we will simply collect all of the targets
  // 
  // For now, only sim 1 trial per condition per person
  array[K, n_trials_to_sim, n_targets] int Q; 
 
  //////////////////////////////////////////////////////////////////////////////
  // first, step through data and compare model selections to human participants
  {
    // some counters and index variables, etc.
    vector[n_targets] weights;  // class weight for teach target
    int t = 0; // counter for trial number
    int kk; // which condition are we in?

    //////////////////////////////////////////////////
    // // step through data row by row and define LLH
    //////////////////////////////////////////////////  
   for (ii in 1:N) {

      t = trial[ii];
      kk = X[t];
   
      // set the weight of each target to be its class weight
      weights = (u_a[kk, Z[ii]]) * to_vector(item_class[t]);

      // multiply weights by stick/switch preference
      weights = inv_logit(weights) .* inv_logit(u_stick[kk, Z[ii]] * S[ii]); 

      // compute spatial weights
      weights = weights .* compute_spatial_weights(found_order[ii], n_targets, 
        u_delta[kk, Z[ii]], u_psi[kk, Z[ii]],
        delta[ii], psi[ii]);
          
      // remove already-selected items, and standarise to sum = 1 
      weights = standarise_weights(weights, n_targets, remaining_items[ii]);   

      P[ii] = categorical_rng(weights);
      W[ii] = weights[Y[ii]];
       
    }
  }
  
}