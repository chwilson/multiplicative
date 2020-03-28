
data {
  int<lower=0> N;
  int<lower=0> K; 
  vector[N] cases; 
  vector[N] day; 
  int[N] county;
  int[N] state; 
}

parameters {
  real B[2];
  real B1_state[K]; 
  real B1_county[KK];
  real B2_state[K];
  real B2_county[KK]; 
  real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  B[1] ~ normal(0,10);
  B[2] ~ exponential(5);
  
  for(i in 1:N)
  log(cases[i]) ~ normal(B[1]+B1_state[state[i]] + B1_county[county[i]] + 
  (B[2] + B2_state[state[i]] + B2_county[county[i]])*day[i], sigma);
}

