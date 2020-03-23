
data {
  int<lower=0> N;
  vector[N] deaths;
  real v_mean;
  real v_scale;
  real R0_mean;
  real R0_scale; 
  real psi_mean;
  real psi_scale; 
  real sigmaI_scale;
  real sigmaD_scale; 
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> v;
  real<lower=0> m;
  real<lower=0> R0; 
  real<lower=0> sigma_I;
  real<lower=0> sigma_D; 
}

transformed parameters{
  real<lower=0> psi = m/v; 
}

model {
  
  
  psi ~ student_t(3,psi_mean,psi_scale);
  v ~ student_t(3,v_mean,v_scale);
  R0 ~ student_t(3,R0_mean,R0_scale);
  I[1] ~ normal();
  
  for(i in 1:N){
  I[i] ~ normal(I[i-1] + (v*(R0-1)-m)*I[i-1],sigma_I);
  Deaths[i] ~ normal(Deaths[i-1] + k*I[i-1], sigma_D);
  }
}

