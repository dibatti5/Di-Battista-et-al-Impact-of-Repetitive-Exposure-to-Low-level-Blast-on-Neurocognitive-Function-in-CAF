functions {
vector derived_cut_points(vector p) {
    int K = num_elements(p);
    vector[K - 1] c;

    real cum_sum = 0;
    for (k in 1:(K - 1)) {
      cum_sum += p[k];
      c[k] = logit(cum_sum);
    }

    return c;
  }
}

data {
    int<lower=0> n;                  // Number of observations
    int<lower=2> G;
    int<lower=1, upper=G> group[n]; 
    int outcome[n];                 // Outcome variable
    int<lower=1> K;                 // Number of categories
}

parameters {
    simplex[K] p; // categoy probabilities
    vector[G-1] gamma_exp; // group index parameter
}

transformed parameters {
  ordered[K-1] cut_points = derived_cut_points(p);
  vector[G] gamma;
  gamma[1] = 0;
  for (g in 2:G)
    gamma[g] = gamma_exp[g - 1];
  
}

model {
  p ~ dirichlet(rep_vector(1, K));
  gamma_exp ~ normal(0, 2);

outcome ~ ordered_logistic(gamma[group], cut_points);
}

generated quantities {
  int y_pred[n];
  for (i in 1:n) {
    y_pred[i] = ordered_logistic_rng(gamma[group[i]], cut_points);
  }
}


