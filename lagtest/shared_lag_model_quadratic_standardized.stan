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

  vector[R] alpha_region;

  real<lower=0> sigma;

  real<lower=0> sigma_lag;

}


transformed parameters {

  vector[K] weight;

  weight = softmax(theta);

}


model {

  // ==================================================
  // Lag smoothing prior
  // ==================================================

  theta[1] ~ normal(0, 1);

  for (k in 2:K) {

    theta[k] - theta[k - 1]
      ~ normal(0, sigma_lag);

  }

  sigma_lag ~ exponential(5);


  // ==================================================
  // Other priors
  // ==================================================

  gamma1 ~ normal(0, 1);

  // Because z is standardized, gamma2 can also use N(0,1)
  gamma2 ~ normal(0, 1);

  alpha_region ~ normal(0, 0.5);

  sigma ~ exponential(1);


  // ==================================================
  // Likelihood
  // ==================================================

  for (n in 1:N) {

    real exposure;

    real z;

    real mu;


    // Z_i(t) = weighted average of lagged temperatures
    exposure =
      dot_product(
        X[n],
        weight
      );


    // Center and standardize latent exposure
    z =
      (exposure - z_center)
      /
      z_scale;


    // Quadratic exposure-response relationship
    mu =
      alpha_region[region[n]]
      +
      gamma1 * z
      +
      gamma2 * square(z);


    y[n] ~ normal(
      mu,
      sigma
    );

  }

}


generated quantities {

  vector[N] log_lik;

  vector[N] mu;

  vector[N] exposure_out;

  vector[N] z_out;

  vector[N] linear_effect;

  vector[N] quadratic_effect;

  vector[N] total_exposure_effect;


  for (n in 1:N) {

    exposure_out[n] =
      dot_product(
        X[n],
        weight
      );


    z_out[n] =
      (
        exposure_out[n]
        -
        z_center
      )
      /
      z_scale;


    linear_effect[n] =
      gamma1 * z_out[n];


    quadratic_effect[n] =
      gamma2 * square(
        z_out[n]
      );


    total_exposure_effect[n] =
      linear_effect[n]
      +
      quadratic_effect[n];


    mu[n] =
      alpha_region[region[n]]
      +
      total_exposure_effect[n];


    log_lik[n] =
      normal_lpdf(
        y[n] |
        mu[n],
        sigma
      );

  }

}
