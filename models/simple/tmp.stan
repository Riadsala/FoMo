
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
  array[N] real log_lik;

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

      // compute spatial weights
      weights = exp(log(weights) +  compute_spatial_weights(found_order[ii], n_targets, 
        rho_delta[kk], rho_psi[kk],
        delta[ii], psi[ii]));
          
      // remove already-selected items, and standarise to sum = 1 
      weights = standarise_weights(weights, n_targets, remaining_items[ii]);   

      P[ii] = categorical_rng(weights);
      log_lik[ii] = weights[Y[ii]];
       
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
          weights = exp(log(weights) + compute_spatial_weights(found_order[jj], n_targets, 
            rho_delta[k], rho_psi[k],
            delta_j, psi_j));
                
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
