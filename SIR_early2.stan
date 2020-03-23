
data {
  int<lower=0> N;
  vector[N] deaths;
  real v_mean;
  real v_scale;
  real R0_mean;
  real R0_scale; 
  real psi_mean;
  real psi_scale; 
  real sigmaD_scale; 
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0> v;
  real<lower=0> m;
  real<lower=0> R0; 
  real<lower=0> I_0; 
  real<lower=0> sigma_D; 
}

transformed parameters{
  real<lower=0> psi = m/v; 
}

model {
  
  sigma_D ~ normal(0, sigmaD_scale); 
  
  I_0 ~ student_t(3,Deaths[1]*750,Deaths[1]*250);
  psi ~ student_t(3,psi_mean,psi_scale);
  v ~ student_t(3,v_mean,v_scale);
  R0 ~ student_t(3,R0_mean,R0_scale);
  D_0 ~ normal(Deaths[1],sigma_D); 
  
  for(i in 2:N){
  Deaths[i] ~ normal(D_0 + I_0*(m/(v*(R_0-1)-m))*(exp((v*(R_0-1)-m)*Day[i])-1), sigma_D);
  }
}
