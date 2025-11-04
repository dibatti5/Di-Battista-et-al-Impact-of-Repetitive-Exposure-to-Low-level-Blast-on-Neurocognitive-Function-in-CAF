//
data {
  int<lower = 1> n; // number of participants
  int<lower = 1> G;
  int<upper = G> group[n]; //group variable
  vector [n] outcome;
  int ord_pred[n];
  int K; // levels of the ordered predictor
  vector [n]age;

}

parameters {
  real<lower = 0> nu;
  vector [G] g;
  real op; 
  real a;
  real<lower = 0> sig;
  simplex[K-1] delta;

}

transformed parameters{
  vector[K] delta_j;
  delta_j = append_row(0, delta);
  
}

model {
  // priors
 g ~ normal( 0 , 0.5 );
 op ~ normal( 0 , 0.5 );
 a ~ normal( 0 , 0.5 );
 nu ~ gamma( 2 , 0.1 );
 sig ~ exponential( 1 );
 delta ~ dirichlet(rep_vector(2, K-1));
 

  //likelihood
  //for (i in 1:n) {
   //outcome[i] ~ student_t(nu, g[group[i]] + a * age[i] + 
               // op * sum(delta_j[1:ord_pred[i]]), sig);
 // }

}

generated quantities{
    vector[n] log_lik;
     vector[n] mu;
    for ( i in 1:n ) {
        mu[i] = g[group[i]] + a * age[i] + 
                op * sum(delta_j[1:ord_pred[i]]);
    }
    for ( i in 1:n ) log_lik[i] = student_t_lpdf( outcome[i] | nu, mu[i] , sig );
}
