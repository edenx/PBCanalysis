// BYM3
data {
  int<lower=0> N;
  int<lower=0> y[N];            // response
  
  vector<lower=0>[N] E;         // exposure
  // matrix[N,7] x;                    // predictor
  // int<lower=0> l;               // number of Fourier features
  matrix[N,N]  L;               // Cholesky decomposition of covariance matrix
  
}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;                 // intercept
  // vector[7] beta1;         // slope

  real<lower=0> sigma;        // overall standard deviation
  real<lower=0, upper=1> rho; // proportion unstructured vs. spatially structured variance

  vector[N] theta;            // unstructured effect
  vector[N] phi;              // spatial effects
}
transformed parameters {
  vector[N] convolved_re;
  // variance of each component should be approximately equal to 1
  convolved_re =  sqrt(1 - rho) * theta + sqrt(rho) * L * phi;
}
model {
// y ~ poisson_log(log_E + beta0 + x * beta1 + convolved_re * sigma);
  y ~ poisson_log(log_E + beta0 + convolved_re * sigma);
  
  beta0 ~ normal(0, 5.0);
  // beta1 ~ normal(0.0, 1.5);
  
  phi ~ std_normal();
  theta ~ normal(0.0, 1.0);
  sigma ~ normal(0, 5);
  rho ~ beta(0.5, 0.5);
}
generated quantities {
  real log_precision = -2.0 * log(sigma);
  real logit_rho = log(rho / (1.0 - rho));
  // vector[N] eta = log_E + beta0 + x * beta1 + convolved_re * sigma;
  vector[N] eta = log_E + beta0 + convolved_re * sigma;
  vector[N] mu = exp(eta);
}

