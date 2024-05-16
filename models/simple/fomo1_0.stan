/* ####### FoMo (Foraging Model 1.0)  #######

This is a new implementation of the model described 
in Clarke et al (2022), Comp Bio.

This version of the model is single level (unclustered data)
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
  array[N] vector<lower = 0>[n_targets] psi; // direction measures (relative)
  array[N] vector[n_targets] phi; // direction measures (absolute)

  array[n_trials] int <lower = 1, upper = K> X; // trial features (ie, which condition are we in)
  array[n_trials, n_targets] int <lower = -1, upper = 1> item_class; // target class, one row per trial
  array[N] vector<lower = -1, upper = 1>[n_targets] S; // stick/switch (does this targ match prev targ) 
  array[N] int<lower = 1, upper = n_trials> trial; // what trial are we on? 

  real prior_sd_b_a; // param for class weight prior
  real prior_sd_b_stick; // prior for sd for b_stick
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
  array[K] real b_a; // weights for class A compared to B  
  array[K] real b_stick; // stick-switch rates 
  array[K] real rho_delta; // distance tuning
  array[K] real rho_psi; // direction tuning 
}

transformed parameters {

}

model {

  // some counters and index variables, etc.
  vector[n_targets] weights;  // class weight for teach target
  vector[n_targets] m; // does this target match the previous target?
  vector[n_targets] prox_weights;
  vector[n_targets] reldir_weights;

  int trl = 0; // counter for trial number
  int kk; // condition (block) index

  /////////////////////////////////////////////////////
  // Define Priors
  ////////////////////////////////////////////////////

  //-----priors intial item selection distributions---

  // priors for fixed effects
  for (ii in 1:K) {
    target += normal_lpdf(b_a[ii]      | 0, prior_sd_b_a);
    target += normal_lpdf(b_stick[ii]      | 0, prior_sd_b_stick);
    target += normal_lpdf(rho_delta[ii] | prior_mu_rho_delta, prior_sd_rho_delta);
    target += normal_lpdf(rho_psi[ii] | prior_mu_rho_psi, prior_sd_rho_psi);
  }

  //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  for (ii in 1:N) {

    // check if we are at the start of a new trial
    // if we are, initialise a load of things
    if (found_order[ii] == 1) {

      trl = trl + 1; // update trial counter           
      kk = X[trl]; // get conditions of current target/trial 

    }

    // new trial, so update the class weights to take random effects into account
    // set the weight of each target to be its class weight
     
    weights = b_a[kk] * to_vector(item_class[trl]);

    // apply spatial weighting
    prox_weights   = compute_prox_weights(found_order[ii], n_targets, 
                                 rho_delta[kk], delta[ii]);
    reldir_weights = compute_reldir_weights(found_order[ii], n_targets, 
                                 rho_psi[kk], psi[ii]);

    if (found_order[ii] == 1) {
      weights = inv_logit(weights);
      } else {
      // check which targets match the previously selected target
      // this is precomputed in S[ii]
      weights = inv_logit(weights) .* inv_logit(b_stick[kk] * S[ii]); 
    }

    weights = weights .* prox_weights .* reldir_weights;

    // remove already-selected items, and standarise to sum = 1
    weights = standarise_weights(weights, n_targets, remaining_items[ii]);

    // likelihood! 
    target += categorical_lpmf(Y[ii] | weights);
    
  }
}

generated quantities {
  // here we  can output our prior distritions
  real prior_b_a = normal_rng(0, prior_sd_b_a);
  real prior_b_stick = normal_rng(0, prior_sd_b_stick);
  real prior_rho_delta = normal_rng(prior_mu_rho_delta, prior_sd_rho_delta);
  real prior_rho_psi = normal_rng(prior_mu_rho_psi, prior_sd_rho_psi);

  /* This code steps through the data item selection by item selection and
  computes
  P: a sampled next target using the model's weights
  W: the weight the model assigned to the target that was selected next

  We will also allow the model to simulate a whole trial start-to-finish so that 
  we can compare run statistics, inter-target selection dynamics etc/
  This will be saved in Q. 
  */
        
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
      weights = b_a[kk] * to_vector(item_class[t]);

      // multiply weights by stick/switch preference
      weights = inv_logit(weights) .* inv_logit(b_stick[kk] * S[ii]); 

      weights = weights .* compute_prox_weights(found_order[ii], n_targets,
         rho_delta[kk], delta[ii]);
          
      // remove already-selected items, and standarise to sum = 1 
      weights = standarise_weights(weights, n_targets, remaining_items[ii]);   

      P[ii] = categorical_rng(weights);
      W[ii] = weights[Y[ii]];
       
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // now allow the model to do a whole trial on its own
  {
    vector[n_targets] remaining_items2;
    vector[n_targets] Sj;
    // some counters and index variables, etc.
    vector[n_targets] weights;  // class weight for teach target

    vector[n_targets] psi_j, phi_j, delta_j;
  
    // for each condition
    for (k in 1:1) {

      //for each trial
      for (t in 1:n_trials_to_sim) {

        // first, set things up!
        remaining_items2 = rep_vector(1, n_targets);

        // simulate a trial!
        for (jj in 1:n_targets) {

          // set the weight of each target to be its class weight
          weights = (b_a[k]) * to_vector(item_class[t]);
           
          // delta_j: distance to previously selected item
          Sj = rep_vector(0, n_targets);
          delta_j = rep_vector(1, n_targets);
          phi_j = rep_vector(1, n_targets);
           
          if (jj > 1) {

            Sj      = compute_matching(item_class[t], n_targets, Q[k, t, ], jj);
            delta_j = compute_prox(item_x[t], item_y[t], n_targets, Q[k, t, ], jj);
            psi_j   = compute_reldir(item_x[t], item_y[t], n_targets, Q[k, t, ], jj); 
              
          }

          // multiply weights by stick/switch preference
          weights = inv_logit(weights) .* inv_logit(b_stick[k] * Sj); 

          // compute spatial weights
          weights = weights .* compute_prox_weights(jj, n_targets,
            rho_delta[k], 
            delta_j);
                
          // remove already-selected items, and standarise to sum = 1 
          weights = standarise_weights(weights, n_targets, remaining_items2);   
     
          Q[k, t, jj] = categorical_rng(weights);

          // update remaining_items2
          remaining_items2[Q[k, t, jj]] = 0;

        }
      }
    }
  }
}
