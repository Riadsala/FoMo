/* 

FoMo V1.3 (multi-level: generated quantities)

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
  int <lower = 1> L; // number of participant levels 
  
  int <lower = 1> n_trials;  // total number of trials (overall)
  int <lower = 1> n_classes; // number of target classes - we assume this is constant over n_trials
  int <lower = 1> n_targets; // total number of targets per trial

  array[N] int <lower = 1> Y; // target IDs - which target was selected here? This is what we predict
  array[n_trials] int <lower = 1, upper = K> X; // trial features (ie, which condition are we in)
  array[n_trials] int <lower = 1, upper = L> Z; // random effect levels

  array[N] int<lower = 1, upper = n_trials> trial; // what trial are we on? 
  array[N] int <lower = 0, upper = n_targets> found_order; // = 1 if starting a new trial, 0 otherwise
  array[n_trials] int first_items;

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

  // hyper parameters
  real<lower = 0> kappa; // kappa = 10? concentration of von Mises
  array[K] real grid_offset; // angular grid offset (ie, should we rotate our coordiante acis)
  real d0; // scale parameter to get rho_delta to ~ 1

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

  // theta is a 4D vector containing the mixture weights for our direction model
  array[K, 4] real log_theta; // mixing values for directions

  ///////////////////////////////
  // random effects
  ///////////////////////////////
  // random effect variances: 
  // 4*K as we have four fixed effect parameters x K conditions
  vector<lower=0>[4*K] sigma_u;
  // 4*K as we have four directions x K conditions
  vector<lower=0>[4*K] sigma_w;
  // declare L_u to be the Choleski factor of a correlation matrix
  cholesky_factor_corr[4*K] L_u;
  // random effect matrix
  matrix[4*K,L] z_u; // for main params
  matrix[4*K, L] z_w; // for directional params
}

transformed parameters {

  /* 
  combine fixed and random effects
  we do this here so that the code in the model{} block is easier to read
  */

  // this transform random effects so that they have the correlation
  // matrix specified by the correlation matrix above
  matrix[4*K, L] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u;

  // create empty arrays for everything
  array[K] vector[L] u_a, u_stick, u_delta, u_psi;
  // calculate
  for (kk in 1:K) {
    u_a[kk]     = to_vector(b_a[kk]       + u[4*(kk-1)+1]);
    u_stick[kk] = to_vector(b_s[kk]   + u[4*(kk-1)+2]);
    u_delta[kk] = to_vector(rho_delta[kk] + u[4*(kk-1)+3]);
    u_psi[kk]   = to_vector(rho_psi[kk]   + u[4*(kk-1)+4]);
  }

  // now work out thet u_log_theta
  array[K, L] vector[4] u_log_theta;

  for (kk in 1:K) {
    for (l in 1:L) {

      u_log_theta[kk, l, 1] = (log_theta[kk, 1] + z_w[4*(kk-1)+1, l] .* sigma_w[4*(kk-1)+1]);
      u_log_theta[kk, l, 2] = (log_theta[kk, 2] + z_w[4*(kk-1)+2, l] .* sigma_w[4*(kk-1)+2]);
      u_log_theta[kk, l, 3] = (log_theta[kk, 3] + z_w[4*(kk-1)+3, l] .* sigma_w[4*(kk-1)+3]);
      u_log_theta[kk, l, 4] = (log_theta[kk, 4] + z_w[4*(kk-1)+4, l] .* sigma_w[4*(kk-1)+4]);

    }  
  }
}

model {

  
}

generated quantities {

  //////////////////////////////////////////////////////////////////////////////
  /* 
  i) Step through data item-selection at a time and compute which item the 
  model selects next. This allows us to compute accuracy

  P: which item do we think the human will select next?
  log_lik: what is the likelihood of selecting the next item ? 
  */
  //////////////////////////////////////////////////////////////////////////////

  array[N] int P;
  array[N] real log_lik;

  {
    
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

      weights = compute_weights_v13(
        u_a[x, z], u_stick[x, z], u_delta[x, z], u_psi[x, z], u_log_theta[x, z], kappa,
        to_vector(item_class[t]), S[ii], delta[ii], psi[ii], phi[ii],
        found_order[ii], n_targets, remaining_items[ii],
        grid_offset[x]); 

      P[ii] = categorical_rng(exp(weights));
      log_lik[ii] = weights[Y[ii]];

    }
  }

  //////////////////////////////////////////////////////////////////////////////
  /* 
  ii) Simulate each trial, start to finish. This allows us to compuare the model
  to human participants in terms of run statistics and selection vectors

  - have not implemented starting rule: select item at random
  - we do not have a stopping rule yet: so we will simply collect all of the targets
  */
  //////////////////////////////////////////////////////////////////////////////

  array[n_trials, n_targets] int Q; 

  {

    // create some variables
    vector[n_targets] weights;
    int z, x; // trial, person, and condition
    /* we need new remaining_items and features (S, psi, delta) as
    these features will update dynamically as we carry out a new 
    simulated foraging trial */
    vector[n_targets] remaining_items_q;
    vector[n_targets] S_q, psi_q, phi_q, delta_q;

    //for each trial
    for (t in 1:n_trials) {

      z = Z[t];
      x = X[t];

      // first, set up new trial_ all the items are remaining!
      remaining_items_q = rep_vector(1, n_targets);
   
      // simulate a trial!
      for (ii in 1:n_targets) {

        // create empty vectors to store our features in
        S_q     = rep_vector(0, n_targets);
        delta_q = rep_vector(1, n_targets);
        psi_q   = rep_vector(1, n_targets);
        phi_q   = rep_vector(1, n_targets);
          
        // if we're not on the first item.... calculate feature vectors  
        if (ii > 1) 
        {
          S_q     = compute_matching(item_class[t], n_targets, Q[t, ], ii);
          delta_q = d0 * compute_prox(item_x[t], item_y[t], n_targets, Q[t, ], ii);
          psi_q   = compute_reldir(item_x[t], item_y[t], n_targets, Q[t, ], ii);
          phi_q   = compute_absdir(item_x[t], item_y[t], n_targets, Q[t, ], ii); 
        
           weights = compute_weights_v13(
          u_a[x, z], u_stick[x, z], u_delta[x, z], u_psi[x, z], u_log_theta[x, z], kappa,
          to_vector(item_class[t]), S_q, delta_q, psi_q, phi_q,
          found_order[ii], n_targets, remaining_items_q,
          grid_offset[x]); 

          // sample an item to select
          Q[t, ii] = categorical_rng(exp(weights));

        } else {

          // first item same as it ever was
          Q[t, ii] = first_items[t];
        }

       

       

        // update remaining_items_q
        remaining_items_q[Q[t, ii]] = 0;
 
      }  
    }
  } 
}
