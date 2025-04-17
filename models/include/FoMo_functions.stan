vector compute_weights_v13(

  // compute weights for FoMo v1.3
  // this is  v1.2 + absolute (up/down/left/right) phi weights

  real b_a, real b_s, real rho_delta, real rho_psi, vector log_theta, real kappa,
  vector item_class, vector match_prev_item, vector delta, vector psi, vector phi, 
  int n, int n_targets, vector remaining_items, 
  real os) {

  // set the weight of each target to be its class weight
  vector[n_targets] weights = log_inv_logit(b_a * to_vector(item_class));

  //  weights by stick/switch preference
  weights += log_inv_logit(b_s * match_prev_item); 

  // apply proximity weights
  weights += -rho_delta * delta - rho_psi * psi;

  // apply abs. dir. weights
  weights += compute_absdir_weights_fixed_kappa4(n, n_targets, 
                                 log_theta, kappa, phi, os);

  // remove already-selected items, and standarise to sum = 1 
  weights += log(remaining_items);
  weights -= log_sum_exp(weights);

  return(weights);

}

vector compute_weights_v12(

  // compute weights for FoMo v1.1
  // this is the same as v1.0 but ignores psi

  real b_a, real b_s, real rho_delta,
  vector item_class, vector match_prev_item, vector delta,
  int n, int n_targets, vector remaining_items) {

  // set the weight of each target to be its class weight
  vector[n_targets] weights = log_inv_logit(b_a * to_vector(item_class));

  //  weights by stick/switch preference
  weights += log_inv_logit(b_s * match_prev_item); 

  // apply proximity weights
  weights += -rho_delta * delta;

  // remove already-selected items, and standarise to sum = 1 
  weights += log(remaining_items);
  weights -= log_sum_exp(weights);

  return(weights);

}

vector compute_weights_v10(

  // compute weights for FoMo v1.0
  // this is our baseline Clarke et al (2022) Comp Bio model

  real b_a, real b_s, real rho_delta, real rho_psi,
  vector item_class, vector match_prev_item, vector delta, vector psi,
  int n, int n_targets, vector remaining_items) {

  // set the weight of each target to be its class weight
  vector[n_targets] weights = log_inv_logit(b_a * to_vector(item_class));

  //  weights by stick/switch preference
  weights += log_inv_logit(b_s * match_prev_item); 

  // apply proximity and rel-dir preference
  weights += -rho_delta * delta - rho_psi * psi;

  // remove already-selected items, and standarise to sum = 1 
  weights += log(remaining_items);
  weights -= log_sum_exp(weights);

  return(weights);

}


////////////////////////////////////////////////
// helper functions
////////////////////////////////////////////////

vector compute_prox(vector x, vector y, int n_targets, array [ ] int Q, int jj) {

  vector[n_targets] delta;

  // compute delta
  delta = (x[Q[jj-1]] - x)^2 + (y[Q[jj-1]] - y)^2;
  delta = sqrt(delta);

  return(delta);

}

array [ ] vector scale_all_prox(array [ ] vector delta, array [ ] vector ri, int N, int n_targets) {
  // normalise delta
  array[N] vector[n_targets] delta_n;

  for (row in 1:N) {

    delta_n[row] = scale_prox(delta[row], ri[row], n_targets);

  }

  return(delta_n);

}

vector scale_prox(vector p, vector ri, int n_targets) {

  real min_delta = 1000;
  vector[n_targets] d;

  // first, set delta to 0 for items that we have already found
  d = p .* ri;

  // we want to find the smallest nonzero value...
  // so let's first set all the 0s to some large value
  for (ii in 1:n_targets) {
    if ((d[ii] < min_delta) && (d[ii] > 0)) {
      min_delta = d[ii];
    }
  }

  // normalise
  d = d ./ min_delta;

  return(d);

}

vector compute_absdir(vector x, vector y, int n_targets, array [ ] int Q, int jj) {      

  vector[n_targets] phi = rep_vector(1, n_targets);

  if (jj > 1) {

    phi = atan2((y - y[Q[jj-1]]), (x - x[Q[jj-1]])) * 180 / pi();

  }

  return(phi);

}


vector compute_reldir(vector x, vector y, int n_targets, array [ ] int Q, int jj) {      

  vector[n_targets] psi = rep_vector(1, n_targets);

  if (jj > 2) {

    psi = atan2((y - y[Q[jj-1]]), (x - x[Q[jj-1]])) * 180 / pi();

    psi = psi - atan2((y[Q[jj-1]] - y[Q[jj-2]]), (x[Q[jj-1]] - x[Q[jj-2]])) * 180 / pi();

    for (tt in 1:n_targets) {
      psi[tt] = fmin(abs(mod(psi[tt], 360)), abs(mod(-psi[tt], 360)))/180;
    }
  }

  return(psi);

}

vector compute_matching(array [ ] int item_class, int n_targets, array [ ] int Q, int jj) {

  // S: which items match the previously selected item?
  vector[n_targets] S;


  for (tt in 1:n_targets) {

   if (item_class[tt] == item_class[Q[jj-1]]) {
    S[tt] = 1;
  } else {
    S[tt] = -1;
  }
}

return(S);
}

real mod(real a, real b) {

  return( a - floor(a/b));

}

vector standarise_weights(vector w, int n_targets, vector remaining_items) {

    /* set weights of found items to 0 and divide by the sum of 
    remaining weights so that they sum to 1 */
    vector[n_targets] w_s = w .* remaining_items;  
    w_s = w_s / sum(w_s);
    return(w_s);
  }
  
  vector compute_absdir_weights_fixed_kappa4(int n, int n_targets, 
    vector log_theta, real kappa, vector phi, real os) {

    // this is a version of the function with only 4 directional components, used in model 1.3, for testing
    // in general, use the version with 8

    vector[n_targets] w;
    vector[4] theta = exp(log_theta);

    w = rep_vector(0, n_targets); 

    if (n > 1) {

      w = w 
      .* (theta[1]*exp(kappa * cos(phi            - os))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[2]*exp(kappa * cos(phi - 2*pi()/4) - os)./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[3]*exp(kappa * cos(phi - 4*pi()/4) - os)./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[4]*exp(kappa * cos(phi - 6*pi()/4) - os)./(2*pi()*modified_bessel_first_kind(0, kappa)) 
       +  1);
    }

    return(w);
  }

  vector compute_absdir_weights_fixed_kappa(int n, int n_targets, 
    vector log_theta, real kappa, vector phi) {

    vector[n_targets] w;
    vector[8] theta = exp(log_theta);

    w = rep_vector(0, n_targets); 

    if (n > 1) {

      w = w 
      .* (theta[1]*exp(kappa * cos(phi           ))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[2]*exp(kappa * cos(phi - 1*pi()/4))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[3]*exp(kappa * cos(phi - 2*pi()/4))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[4]*exp(kappa * cos(phi - 3*pi()/4))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[5]*exp(kappa * cos(phi - 4*pi()/4))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[6]*exp(kappa * cos(phi - 5*pi()/4))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[7]*exp(kappa * cos(phi - 6*pi()/4))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  theta[8]*exp(kappa * cos(phi - 7*pi()/4))./(2*pi()*modified_bessel_first_kind(0, kappa))
       +  1);
    }

    return(w);
  }

  vector compute_prox_weights(int n, int n_targets, real rho_delta, vector delta) {

    vector[n_targets] w = -rho_delta * delta;
    return(w);
  }

  vector compute_reldir_weights(int n, int n_targets, real rho_psi, vector psi) {

    vector[n_targets] w = - rho_psi * psi;
    return(w);
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
        remaining_items[n][Y[n-1]] = 0.0000;

      }
    }

    return(remaining_items);
  }
