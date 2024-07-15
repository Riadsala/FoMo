/* ####### FoMo (Foraging Model 0)  #######

Simple test version

*/

functions {

vector compute_prox_weights(int n, int n_targets, real rho_delta, vector delta) {

  vector[n_targets] w;

  w = rep_vector(1, n_targets); 
  
  // now start computing the weights
  if (n > 1) {
    // for the second selected target, weight by distance from the first
    w = w .* exp(-rho_delta * delta);
    
  }

  return(w);
}


  vector standarise_weights(vector w, int n_targets, vector remaining_items) {

    /* set weights of found items to 0 and divide by the sum of 
    remaining weights so that they sum to 1 */
    vector[n_targets] w_s = w .* remaining_items;  
    w_s = w_s / sum(w_s);
    return(w_s);
  }
  
  array[ ] vector calc_remaining_items(int N, int n_targets, array[ ] int Y, array[ ] int found_order) {

  array[N] vector[n_targets] remaining_items;

  for (n in 1:N) {
    // check if we are at the start of a new trial
    // if we are, initialise a load of things
    if (found_order[n] == 1) {
             
      // as we're at the start of a new trial, reset the remaining_items tracker
      remaining_items[n] = rep_vector(1, n_targets);

    } else {

      remaining_items[n] = remaining_items[n-1];
      remaining_items[n][Y[n-1]] = 0;
      
    }
  }

  return(remaining_items);
}

  vector compute_spatial_weights(int n, int n_targets, 
    real rho_delta, vector delta) {

    // computes spatial weights
    // for FoMo0, this includes proximity only
    vector[n_targets] prox_weights;

    // apply spatial weighting
    prox_weights  = compute_prox_weights(n, n_targets, 
                                 rho_delta, delta);
   
    // return the dot product of the weights
    return(prox_weights);

  }
}

data {
  int <lower = 1> N; // total number of selected targets over the whole experiment
  int <lower = 1> K; // number of experimental conditions  

  int <lower = 1> n_trials;  // total number of trials (overall)
  int <lower = 1> n_targets; // total number of targets per trial
  array[N] int <lower = 0, upper = n_targets> found_order; // = 1 is starting a new trial, 0 otherwise

  array[N] int <lower = 1> Y; // target IDs - which target was selected here? This is what we predict
  array[N] vector<lower = 0>[n_targets] delta; // distance measures

  array[N] int<lower = 1, upper = n_trials> trial; // what trial are we on? 
 
  real prior_mu_rho_delta;
  real prior_sd_rho_delta;
}

transformed data{
  // compute remaining items
  array[N] vector[n_targets] remaining_items;
  remaining_items = calc_remaining_items(N, n_targets, Y, found_order);
  
}

parameters {
  
  // These are all the parameters we want to fit to the data
  real rho_delta; // distance tuning

}

model {

  // vector to store weights for each target
  vector[n_targets] weights;

  /////////////////////////////////////////////////////
  // Define Priors
  ////////////////////////////////////////////////////
  target += normal_lpdf(rho_delta | prior_mu_rho_delta, prior_sd_rho_delta);
    
   //////////////////////////////////////////////////
  // // step through data row by row and define LLH
  //////////////////////////////////////////////////  
  for (ii in 1:N) {

    // new trial, so update the class weights to take random effects into account
    // set the weight of each target to be its class weight

    // apply spatial weighting
    weights = compute_spatial_weights(found_order[ii], n_targets, 
      rho_delta, delta[ii]);

    // remove already-selected items, and standarise to sum = 1
    weights = standarise_weights(weights, n_targets, remaining_items[ii]);

    // likelihood! 
    target += categorical_lpmf(Y[ii] | weights);
    
  }
}

generated quantities {
  // here we  can output our prior distritions
  real prior_rho_delta = normal_rng(prior_mu_rho_delta, prior_sd_rho_delta);
 

  /* This code steps through the data item selection by item selection and
  computes
  P: a sampled next target using the model's weights
  W: the weight the model assigned to the target that was selected next

  We will also allow the model to simulate a whole trial start-to-finish so that 
  we can compare run statistics, inter-target selection dynamics etc/
  This will be saved in Q. 
  */
        
  array[N] real log_lik;

  // for trial level predictions, we have to remember that we do not have a stopping rule yet
  // so we will simply collect all of the targets
 
  //////////////////////////////////////////////////////////////////////////////
  // first, step through data and compare model selections to human participants

    // some counters and index variables, etc.
    vector[n_targets] weights;  // class weight for teach target
 
    //////////////////////////////////////////////////
    // // step through data row by row and define LLH
    //////////////////////////////////////////////////  
   for (ii in 1:N) {

      // compute spatial weights
      weights = compute_spatial_weights(found_order[ii], n_targets, 
        rho_delta, delta[ii]);
        
      // remove already-selected items, and standarise to sum = 1 
      weights = standarise_weights(weights, n_targets, remaining_items[ii]);   

      log_lik[ii] = weights[Y[ii]];
       
    }
  }

