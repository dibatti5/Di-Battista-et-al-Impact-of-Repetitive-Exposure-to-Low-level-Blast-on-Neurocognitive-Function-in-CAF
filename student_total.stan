//
data {
  int<lower = 0> n; // number of participants
  int<lower = 1> G;
  int<lower = 1, upper = G> group[n]; //group variable
  vector [n]outcome;
}

parameters {
  real<lower = 0> nu;
  vector [G] g;
  real<lower = 0> sig;
}

model {
  // priors
 g ~ normal( 0 , 0.4 );
 nu ~ gamma( 2 ,0.1 );
 sig ~ exponential( 1 );

  //likelihood
  for (i in 1:n) {
    outcome[i] ~ student_t(nu, g[group[i]], sig);
  }
}

generated quantities{
    vector[n] log_lik;
     vector[n] mu;
    for ( i in 1:n ) {
        mu[i] = g[group[i]];
    }
    for ( i in 1:n ) log_lik[i] = student_t_lpdf( outcome[i] | nu, mu[i] , sig );
}
