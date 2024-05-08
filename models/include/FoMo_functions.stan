vector compute_prox(vector x, vector y, int n_targets, array [ ] int Q, int jj) {

  vector[n_targets] delta;

  // compute delta
  delta = (x[Q[jj-1]] - x)^2 + (y[Q[jj-1]] - y)^2;
  delta = sqrt(delta);

  return(delta);

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

  vector compute_spatial_weights(int n, int n_targets, 
                                 real rho_delta, real rho_psi,
                                 vector delta, vector psi, vector phi,
                                 vector item_x, vector item_y) {

    vector[n_targets] w;

    w = rep_vector(1, n_targets); 
     // now start computing the weights
    if (n == 1) {

      // calculate inital selection weights based on spatial location

    } else {

      if (n == 2) {
        // for the second selected target, weight by distance from the first
        w = w .* exp(-rho_delta * delta);
      } else {
        // for all later targets, also weight by direciton
        w = w .* exp(-rho_delta * delta - rho_psi * psi);
      }
    }
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
      remaining_items[n][Y[n-1]] = 0;
      
    }
  }

  return(remaining_items);
}
