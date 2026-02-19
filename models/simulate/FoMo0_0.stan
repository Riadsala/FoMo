/* 

FoMo v0.0 (closest item baseline)

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

  ///////////////////////////////
  // random effects
  ///////////////////////////////
  // random effect variances: 
  // 4*K as we have four fixed effect parameters x K conditions
  vector<lower=0>[4*K] sigma_u;
  // declare L_u to be the Choleski factor of a correlation matrix
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


      weights = (inv(delta[ii]) .* remaining_items[ii]);

      real max_weights = max(weights);

      int p = 0;

      for (j in 1:n_targets) {
        if (weights[j] == max_weights) {
          p = j;
        }

        P[ii] = p;
      }


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

}
