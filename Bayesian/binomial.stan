data {
  int<lower=1> T;
  int k[T];
  int n[T];
}
transformed data {
}
parameters {
  real<lower=0, upper=1> theta;
}
transformed parameters {
}
model {
  theta ~ beta(2, 2);
  for (t in 1:T) {
    k[t] ~ binomial(n[t], theta);
  }
}
generated quantities{
}
