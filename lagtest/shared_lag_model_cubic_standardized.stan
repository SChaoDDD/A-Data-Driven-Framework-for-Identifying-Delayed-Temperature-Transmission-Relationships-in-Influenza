data {

  int<lower=1> N;
  vector[N] y;

  int<lower=1> K;
  matrix[N, K] X;

  int<lower=1> R;
  array[N] int<lower=1, upper=R> region;

  // Centering and scaling constants for latent exposure Z
  real z_center;
  real<lower=0> z_scale;

}


parameters {

  vector[K] theta;

  real gamma1;
  real gamma2;
  real gamma3;

  vector[R] alpha_region;

  real<lower=0> sigma;
  real<lower=0> sigma_lag;

}


transformed parameters {

  vector[K] weight;

  weight = softmax(theta);

}


model {

  // ======================
  // lag smoothing prior
  // ======================

  theta[1] ~ normal(0, 1);

  for (k in 2:K) {
    theta[k] - theta[k - 1]
      ~ normal(0, sigma_lag);
  }

  sigma_lag ~ exponential(5);


  // ======================
  // other priors
  // ======================

  gamma1 ~ normal(0, 1);
  gamma2 ~ normal(0, 1);
  gamma3 ~ normal(0, 1);

  alpha_region ~ normal(0, 0.5);

  sigma ~ exponential(1);


  // ======================
  // likelihood
  // ======================

  for (n in 1:N) {

    real exposure;
    real z;
    real mu;

    exposure =
      dot_product(
        X[n],
        weight
      );

    z =
      (exposure - z_center)
      /
      z_scale;

    mu =
      alpha_region[region[n]]
      +
      gamma1 * z
      +
      gamma2 * square(z)
      +
      gamma3 * z * z * z;

    y[n] ~ normal(
      mu,
      sigma
    );

  }

}
