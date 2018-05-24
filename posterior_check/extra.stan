// model with alpha and beta
data {
  int<lower=1> N;
  int<lower=1> T;               
  int<lower=1,upper=T> Tsubj[N];                 
  int<lower=-1,upper=2> choice[N,T];     
  real outcome[N,T];  // no lower and upper bounds   
}

transformed data {
  vector[2] initV;  // initial values for EV
  initV = rep_vector(0.0, 2);
}

parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters  
  vector[3] mu_p;  
  vector<lower=0>[3] sigma;
    
  // Subject-level raw parameters (for Matt trick)
  vector[N] alpha_pr;    // learning rate [0, 1]
  vector[N] beta_pr;  // inverse temperature [0, 5]
}

transformed parameters {
  // subject-level parameters
  vector<lower=0,upper=1>[N] alpha_pos;
  vector<lower=0,upper=1>[N] alpha_neg;
  vector<lower=0,upper=5>[N] beta;
  
  for (i in 1:N) {
    alpha_pos[i] = Phi_approx( mu_p[1] + sigma[1] * alpha_pr[i] );
    alpha_neg[i] = Phi_approx( mu_p[2] + sigma[2] * alpha_pr[i] );
    beta[i] = Phi_approx( mu_p[3] + sigma[3] * beta_pr[i] ) * 5;
  }
}

model {
  // Hyperparameters
  mu_p  ~ normal(0, 1); 
  sigma ~ normal(0, 1);  
  
  // individual parameters
  alpha_pr ~ normal(0,1);
  beta_pr ~ normal(0,1);
  
  // subject loop and trial loop
  for (i in 1:N) {
    vector[2] ev; // expected value
    real PE;      // prediction error

    ev = initV;

    for (t in 1:(Tsubj[i])) {        
      // compute action probabilities
      choice[i,t] ~ categorical_logit( beta[i]*ev );

      // prediction error 
      PE = outcome[i,t] - ev[choice[i,t]];

      // value updating (learning) 
      if(PE >= 0) {
        ev[choice[i,t]] = ev[choice[i,t]] + alpha_pos[i] * PE; 
      } else {
        ev[choice[i,t]] = ev[choice[i,t]] + alpha_neg[i] * PE; 
      }
    }
  }
}
generated quantities {
  real<lower=0,upper=1> mu_alpha_pos;
  real<lower=0,upper=1> mu_alpha_neg;
  real<lower=0,upper=5> mu_beta;
  
  mu_alpha_pos = Phi_approx(mu_p[1]);
  mu_alpha_neg = Phi_approx(mu_p[2]);
  mu_beta = Phi_approx(mu_p[3]) * 5;
}