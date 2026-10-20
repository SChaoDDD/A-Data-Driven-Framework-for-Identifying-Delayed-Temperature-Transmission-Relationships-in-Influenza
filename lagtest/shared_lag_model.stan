data {

  int<lower=1> N;

  vector[N] y;

  int<lower=1> K;

  matrix[N,K] X;

  int<lower=1> R;

  array[N] int region;

}


parameters {

  vector[K] theta;

  real gamma;

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


  // 第一阶lag没有特殊偏好
  theta[1] ~ normal(0,1);


  // 相邻lag平滑
  for(k in 2:K){

    theta[k]-theta[k-1]
      ~ normal(0,sigma_lag);

  }


  sigma_lag ~ exponential(5);



  // ======================
  // other priors
  // ======================


  gamma ~ normal(0,1);


  alpha_region ~ normal(0,0.5);


  sigma ~ exponential(1);



  // ======================
  // likelihood
  // ======================


  for(n in 1:N){

    real exposure;


    exposure =
      dot_product(
        X[n],
        weight
      );


    y[n]
    ~ normal(
      alpha_region[region[n]]
      +
      gamma*exposure,
      sigma
    );

  }

}

generated quantities {

  vector[N] log_lik;


  for(n in 1:N){

    real exposure;


    exposure =
      dot_product(
        X[n],
        weight
      );


    log_lik[n] =
      normal_lpdf(
        y[n] |
        alpha_region[region[n]]
        + gamma * exposure,
        sigma
      );

  }

}
