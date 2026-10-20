# ============================================================
# Bayesian shared lag model
# Centered + standardized quadratic exposure-response
#
# Z_i(t) = sum_l w_l * Temp_i(t-l)
#
# z_i(t) = (Z_i(t) - z_center) / z_scale
#
# y_i(t) ~ Normal(
#   alpha_i + gamma1*z_i(t) + gamma2*z_i(t)^2,
#   sigma^2
# )
#
# Stan model:
#   shared_lag_model_quadratic_standardized.stan
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(cmdstanr)
library(loo)
library(rstan)

rstan_options(auto_write = TRUE)

options(
  mc.cores = parallel::detectCores()
)


# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

base_dir <-
  "C:/Users/XLYC/Desktop/流感代码修改/AAAinfluenzacode/lagtest"


ia_file <-
  file.path(
    base_dir,
    "IAbeta.csv"
  )

il_file <-
  file.path(
    base_dir,
    "ILbeta.csv"
  )

in_file <-
  file.path(
    base_dir,
    "INbeta.csv"
  )

stan_file <-
  file.path(
    base_dir,
    "shared_lag_model_quadratic_standardized.stan"
  )


# ------------------------------------------------------------
# 2. Check files
# ------------------------------------------------------------

if (!file.exists(ia_file)) {
  stop(
    "Cannot find: ",
    ia_file
  )
}

if (!file.exists(il_file)) {
  stop(
    "Cannot find: ",
    il_file
  )
}

if (!file.exists(in_file)) {
  stop(
    "Cannot find: ",
    in_file
  )
}

if (!file.exists(stan_file)) {
  stop(
    "Cannot find Stan file: ",
    stan_file
  )
}


# ------------------------------------------------------------
# 3. Read data
# ------------------------------------------------------------

IA <- read.csv(
  ia_file
)

IL <- read.csv(
  il_file
)

IN <- read.csv(
  in_file
)


# ------------------------------------------------------------
# 4. Region labels
# ------------------------------------------------------------

IA$region <- "IA"

IL$region <- "IL"

IN$region <- "IN"


# ------------------------------------------------------------
# 5. Merge data
# ------------------------------------------------------------

data <- bind_rows(
  IA,
  IL,
  IN
)


# ------------------------------------------------------------
# 6. Basic checks
# ------------------------------------------------------------

required_columns <-
  c(
    "Day",
    "beta",
    "temp",
    "region"
  )


missing_columns <-
  setdiff(
    required_columns,
    names(data)
  )


if (length(missing_columns) > 0) {

  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )

}


if (
  any(
    data$beta <= 0,
    na.rm = TRUE
  )
) {

  stop(
    "beta must be > 0 because y = log(beta / beta_mean)."
  )

}


# ------------------------------------------------------------
# 7. Construct response variable
#
# y_i(t) = log(beta_i(t) / regional mean beta)
# ------------------------------------------------------------

data <-
  data %>%
  arrange(
    region,
    Day
  ) %>%
  group_by(
    region
  ) %>%
  mutate(

    beta_mean =
      mean(
        beta,
        na.rm = TRUE
      ),

    y =
      log(
        beta /
          beta_mean
      )

  ) %>%
  ungroup()


# ------------------------------------------------------------
# 8. Maximum lag
# ------------------------------------------------------------

L <- 7


lag_names <-
  paste0(
    "temp_lag",
    0:L
  )


# ------------------------------------------------------------
# 9. Create lagged temperature variables
# ------------------------------------------------------------

for (i in 0:L) {

  data[[
    paste0(
      "temp_lag",
      i
    )
  ]] <-

    ave(

      data$temp,

      data$region,

      FUN =
        function(x) {

          dplyr::lag(
            x,
            n = i
          )

        }

    )

}


# ------------------------------------------------------------
# 10. Remove incomplete lag observations
# ------------------------------------------------------------

model_vars <-
  c(
    "Day",
    "region",
    "beta",
    "temp",
    "y",
    lag_names
  )


keep <-
  complete.cases(
    data[, model_vars]
  )


data_model <-
  data[
    keep,
    ,
    drop = FALSE
  ]


if (nrow(data_model) == 0) {

  stop(
    "No complete observations remain after creating lag variables."
  )

}


cat(
  "Number of observations used:",
  nrow(data_model),
  "\n"
)


# ------------------------------------------------------------
# 11. Construct X
# ------------------------------------------------------------

X <-
  data_model %>%
  select(
    all_of(
      lag_names
    )
  ) %>%
  as.matrix()


storage.mode(X) <-
  "double"


y <-
  as.numeric(
    data_model$y
  )


region_factor <-
  factor(
    data_model$region,
    levels =
      c(
        "IA",
        "IL",
        "IN"
      )
  )


region_id <-
  as.integer(
    region_factor
  )


K <-
  ncol(
    X
  )


R <-
  nlevels(
    region_factor
  )


# ------------------------------------------------------------
# 12. Center and scale the latent exposure Z
#
# Since weight = softmax(theta):
#
#   weight_l > 0
#   sum_l weight_l = 1
#
# Z = X %*% weight is a weighted average of lagged temperatures.
#
# We therefore use the pooled lag-temperature distribution
# to define a fixed center and scale.
# ------------------------------------------------------------

x_pool <-
  as.numeric(
    X
  )


z_center <-
  mean(
    x_pool,
    na.rm = TRUE
  )


z_scale <-
  sd(
    x_pool,
    na.rm = TRUE
  )


if (
  !is.finite(
    z_scale
  )
  ||
  z_scale <= 0
) {

  stop(
    "Temperature SD is zero or invalid."
  )

}


cat(
  "z_center =",
  z_center,
  "\n"
)

cat(
  "z_scale =",
  z_scale,
  "\n"
)


# ------------------------------------------------------------
# 13. Stan data
# ------------------------------------------------------------

stan_data <-
  list(

    N =
      nrow(
        data_model
      ),

    y =
      y,

    K =
      K,

    X =
      X,

    R =
      R,

    region =
      region_id,

    z_center =
      z_center,

    z_scale =
      z_scale

  )


# ------------------------------------------------------------
# 14. Compile Stan model
# ------------------------------------------------------------

model <-
  cmdstan_model(
    stan_file,
    force_recompile = FALSE
  )


# ------------------------------------------------------------
# 15. MCMC sampling
#
# For debugging, you can first use 500 + 500.
# After confirming the model runs well, restore 2000 + 2000.
# ------------------------------------------------------------

fit <-
  model$sample(

    data =
      stan_data,

    seed =
      123,

    chains =
      4,

    parallel_chains =
      4,

    iter_warmup =
      2000,

    iter_sampling =
      2000,

    refresh =
      100

  )



# ------------------------------------------------------------
# 17. PSIS-LOO
# ------------------------------------------------------------

log_lik <-
  fit$draws(
    variables =
      "log_lik"
  )


log_lik_matrix <-
  posterior::as_draws_matrix(
    log_lik
  )


loo_result <-
  loo(
    log_lik_matrix
  )


print(
  loo_result
)


write.csv(

  as.data.frame(
    loo_result$estimates
  ),

  file.path(
    base_dir,
    paste0(
      "loo_quadratic_standardized_L",
      L,
      ".csv"
    )
  ),

  row.names =
    TRUE

)



# ------------------------------------------------------------
# 18. Lag weights
# ------------------------------------------------------------

weight_summary <-
  fit$summary(
    
    variables = "weight",
    
    posterior::default_summary_measures(),
    
    q025 = ~ posterior::quantile2(
      .,
      probs = 0.025
    ),
    
    q50 = ~ posterior::quantile2(
      .,
      probs = 0.50
    ),
    
    q975 = ~ posterior::quantile2(
      .,
      probs = 0.975
    )
    
  )


weight_summary$lag <- 0:L

write.csv(
  
  weight_summary,
  
  file.path(
    base_dir,
    "weight7_2.csv"
    
  ),
  
  row.names =
    FALSE
  
)

