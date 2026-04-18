

data {

  int<lower=1> N;

  int<lower=1> K;

  int<lower=1> G;

  int<lower=1,upper=3> y[N];

  matrix[N, K] X;

  int<lower=1,upper=G> ideology_id[N];

}

parameters {

  vector[K] beta;

  ordered[2] c;

  real mu_alpha;

  real<lower=0> sigma_alpha;

  vector[G] z_alpha;

}

transformed parameters {

  vector[G] alpha;

  alpha = mu_alpha + sigma_alpha * z_alpha;

}

model {

  beta ~ normal(0, 1);

  c ~ normal(0, 2.5);

  mu_alpha ~ normal(0, 1);

  sigma_alpha ~ normal(0, 1);

  z_alpha ~ normal(0, 1);

  for (n in 1:N) {

    y[n] ~ ordered_logistic(alpha[ideology_id[n]] + X[n] * beta, c);

  }

}

generated quantities {

  vector[N] log_lik;

  int<lower=1,upper=3> y_rep[N];

  for (n in 1:N) {

    real eta_n;

    eta_n = alpha[ideology_id[n]] + X[n] * beta;

    log_lik[n] = ordered_logistic_lpmf(y[n] | eta_n, c);

    y_rep[n] = ordered_logistic_rng(eta_n, c);

  }

}


