data {
  int<lower=1> K;
  int<lower=0> N;
  int<lower=0,upper=1> y[N];
  matrix[N, K] x;
}

parameters {
  vector[K] beta;
}

model {
  y ~ bernoulli(Phi(x * beta));
}
