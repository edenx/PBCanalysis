// RFF for given bandwidth
data {
  int<lower=1> n; // number of samples (cells)
  int<lower=1> m; // number of fourier features
  matrix[n,2] x;  // centroid of the computational grid
  vector[n] y;  // count per cell
  vector[n] w;  // normalised population density per cell
  real<lower=0> bw; // bandwidth (?)
  matrix[m,2] omega;  // sampled from kernel corresponded sdf with dim m by d(=2)
  
  // following the recommended prior of betas for regularisation
  real<lower=0> nu; // df of the prior on scale parameter of betas
  real<lower=0> sigma; // scale of the prior ...
}
transformed data {
  matrix[n,m] cosfeatures;
  matrix[n,m] sinfeatures;
  matrix[n,2*m] features;
  
  features = x * omega' * bw;
  
  for(i in 1:n)
    for(j in 1:m) {
      cosfeatures[i,j] = cos(features[i,j]);
      sinfeatures[i,j] = sin(features[i,j]);
    }
  cosfeatures = cosfeatures .* w;
  sinfeatures = sinfeatures .* w;
}
parameters {
  vector[m] alpha;
  vector[m] beta1;
  vector[m] beta2;
  real<lower=0> sigma2;
  real<lower=0> tao;
}
transformed parameters {
  vector[n] fhat;
  fhat = alpha + cosfeatures * beta1 + sinfeatures * beta2;
}
model {
  alpha ~ normal(0, 10);
  beta1 ~ normal(0, tao);
  beta2 ~ normal(0, tao);
  tao ~ student_t(nu, 0, sigma)
  sigma2 ~ normal(0, 1);
  
  y ~ normal(fhat, sigma2);
}
