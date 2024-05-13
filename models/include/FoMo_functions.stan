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
