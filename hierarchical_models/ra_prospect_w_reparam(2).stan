data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  int<lower=-1, upper=1> gamble[N,T];
  real<lower=0> gain[N,T];
  real cert[N,T];  
  real<lower=0> loss[N,T];  // absolute loss amount
}
transformed data {
}
parameters {
  vector[3] mu_p;  
  vector<lower=0>[3] sigma;  
  vector[N] rho_p;  
  vector[N] lambda_p;  
  vector[N] tau_p;  
}
transformed parameters {
  vector<lower=0,upper=2>[N] rho;
  vector<lower=0,upper=5>[N] lambda;
  vector<lower=0>[N] tau;
     
  for (i in 1:N) {
    rho[i]    = Phi_approx( mu_p[1] + sigma[1] * rho_p[i] ) * 2;
    lambda[i] = Phi_approx( mu_p[2] + sigma[2] * lambda_p[i] ) * 5; 
  }
  tau = exp( mu_p[3] + sigma[3] * tau_p );
}
model {
  // ra_prospect: Original model in Soko-Hessner et al 2009 PNAS
  // hyper parameters
  mu_p  ~ normal(0, 1.0); 
  sigma ~ normal(0, 1.0);
  
  // individual parameters w/ Matt trick
  rho_p    ~ normal(0, 1.0);   
  lambda_p ~ normal(0, 1.0);   
  tau_p    ~ normal(0, 1.0);
  
  for (i in 1:N) {
    for (t in 1:Tsubj[i]) {
      real evSafe;    // evSafe, evGamble, pGamble can be a scalar to save memory and increase speed. 
      real evGamble;  // they are left as arrays as an example for RL models.
      real pGamble;
      
      // loss[i,t]=absolute amount of loss (pre-converted in R)
      evSafe   = pow(cert[i,t], rho[i]);    
      evGamble = 0.5 * (pow(gain[i,t], rho[i]) - lambda[i] * pow(loss[i,t], rho[i]) ); 
      pGamble  = inv_logit( tau[i] * (evGamble - evSafe) );
      gamble[i,t] ~ bernoulli( pGamble );
    }
  }
}
