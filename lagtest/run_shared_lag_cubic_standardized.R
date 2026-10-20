# ============================================================
# Bayesian shared lag model
# Centered + standardized cubic exposure-response
#
# Z_i(t) = sum_l w_l * Temp_i(t-l)
#
# z_i(t) = (Z_i(t) - z_center) / z_scale
#
# y_i(t) ~ Normal(
#   alpha_i
#   + gamma1*z_i(t)
#   + gamma2*z_i(t)^2
#   + gamma3*z_i(t)^3,
#   sigma^2
# )
#
# Only posterior weight output is retained.
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

library(dplyr)
library(cmdstanr)
library(posterior)

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
    "shared_lag_model_cubic_standardized.stan"
  )


# ------------------------------------------------------------
# 2. Read data
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
# 3. Region labels
# ------------------------------------------------------------

IA$region <- "IA"
IL$region <- "IL"
IN$region <- "IN"


# ------------------------------------------------------------
# 4. Merge data
# ------------------------------------------------------------

data <- bind_rows(
  IA,
  IL,
  IN
)


# ------------------------------------------------------------
# 5. Response variable
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
# 6. Maximum lag
# ------------------------------------------------------------

L <- 7


lag_names <-
  paste0(
    "temp_lag",
    0:L
  )


# ------------------------------------------------------------
# 7. Create lagged temperature variables
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
# 8. Remove incomplete observations
# ------------------------------------------------------------

model_vars <-
  c(
    "y",
    "region",
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


# ------------------------------------------------------------
# 9. Construct X, y, region
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
    data_model$region
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
# 10. Center and standardize latent exposure Z
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
  !is.finite(z_scale)
  ||
  z_scale <= 0
) {

  stop(
    "Temperature SD is zero or invalid."
  )

}


# ------------------------------------------------------------
# 11. Stan data
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
# 12. Compile Stan model
# ------------------------------------------------------------

model <-
  cmdstan_model(
    stan_file,
    force_recompile = FALSE
  )


# ------------------------------------------------------------
# 13. Sampling
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
# 14. Extract ONLY weight posterior draws
# ------------------------------------------------------------

weight_draws <-
  fit$draws(

    variables =
      "weight",

    format =
      "draws_matrix"

  )


weight_draws <-
  as.matrix(
    weight_draws
  )


# ------------------------------------------------------------
# 15. Weight summary
# ------------------------------------------------------------

weight_summary <-
  data.frame(

    lag =
      0:L,

    variable =
      colnames(
        weight_draws
      ),

    mean =
      apply(
        weight_draws,
        2,
        mean
      ),

    sd =
      apply(
        weight_draws,
        2,
        sd
      ),

    q025 =
      apply(
        weight_draws,
        2,
        quantile,
        probs = 0.025
      ),

    q50 =
      apply(
        weight_draws,
        2,
        quantile,
        probs = 0.50
      ),

    q975 =
      apply(
        weight_draws,
        2,
        quantile,
        probs = 0.975
      )

  )


print(
  weight_summary
)


# ------------------------------------------------------------
# 16. Save ONLY weight output
# ------------------------------------------------------------

write.csv(

  weight_summary,

  file.path(
    base_dir,
      "weight7_3.csv"

  ),

  row.names =
    FALSE

)
