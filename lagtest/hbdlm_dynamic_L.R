
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(cmdstanr)
library(loo)
library(rstan) 
rstan_options(auto_write = TRUE)

options(mc.cores = parallel::detectCores())
IA <- read.csv("C:/Users/XLYC/Desktop/流感代码修改/AAAinfluenzacode/lagtest/IAbeta.csv")
IL <- read.csv("C:/Users/XLYC/Desktop/流感代码修改/AAAinfluenzacode/lagtest/ILbeta.csv")
IN <- read.csv("C:/Users/XLYC/Desktop/流感代码修改/AAAinfluenzacode/lagtest/INbeta.csv")


IA$region <- "IA"
IL$region <- "IL"
IN$region <- "IN"


data <- bind_rows(
  IA,
  IL,
  IN
)


data <- data %>%
  arrange(region, Day)
data <- data %>%
  group_by(region) %>%
  mutate(
    
    beta_mean = mean(beta),
    
    y =
      log(beta / beta_mean)
    
  ) %>%
  ungroup()
# 最大滞后天数，只需要修改这里
L <- 28

# 根据L自动生成temp_lag0 ~ temp_lagL
lag_names <- paste0("temp_lag", 0:L)

data <- data %>%
  group_by(region) %>%
  arrange(Day) %>%
  mutate(
    across(
      .cols = temp,
      .fns = list(),
      .names = "{.col}"
    )
  ) %>%
  ungroup()

# 动态创建lag变量
for (i in 0:L) {
  data[[paste0("temp_lag", i)]] <- ave(
    data$temp,
    data$region,
    FUN = function(x) dplyr::lag(x, i)
  )
}


# 删除没有lag的数据

data <- na.omit(data)
X <- data %>%
  select(all_of(lag_names)) %>%
  as.matrix()
y <- data$y
region_id <- as.integer(
  factor(data$region)
)
stan_data <- list(
  
  N = nrow(data),
  
  y = y,
  
  K = ncol(X),
  
  X = X,
  
  R = length(unique(region_id)),
  
  region = region_id
  
)


model <- cmdstan_model(
  "C:/Users/XLYC/Desktop/流感代码修改/AAAinfluenzacode/lagtest/shared_lag_model.stan"
)


fit <- model$sample(
  
  data=stan_data,
  
  seed=123,
  
  chains=4,
  
  parallel_chains=4,
  
  iter_warmup=2000,
  
  iter_sampling=2000
  
)


# 提取log likelihood

log_lik <- fit$draws(
  variables = "log_lik"
)


# 转换成矩阵

log_lik_matrix <- posterior::as_draws_matrix(log_lik)



# 计算PSIS-LOO

loo_result <- loo(
  log_lik_matrix
)


print(loo_result)
weight21=fit$summary(
  variables="weight"
)
weight21
write.csv(weight21, "weight28_1.csv", row.names = FALSE)
