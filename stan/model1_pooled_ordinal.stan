data {
  int<lower=1> N;
  int<lower=1> K;
  int<lower=1,upper=3> y[N];
  matrix[N, K] X;
}

parameters {
  vector[K] beta;
  ordered[2] c;
}

model {
  beta ~ normal(0, 1);
  c ~ normal(0, 2.5);
  y ~ ordered_logistic(X * beta, c);
}

generated quantities {
  vector[N] log_lik;
  int<lower=1,upper=3> y_rep[N];
  for (n in 1:N) {
    log_lik[n] = ordered_logistic_lpmf(y[n] | (X[n] * beta), c);
    y_rep[n] = ordered_logistic_rng((X[n] * beta), c);
  }
}


