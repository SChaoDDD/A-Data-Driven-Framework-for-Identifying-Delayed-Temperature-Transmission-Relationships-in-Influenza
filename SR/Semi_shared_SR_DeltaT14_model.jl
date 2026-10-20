
# 第一次运行时取消下面几行注释
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("Statistics")
# Pkg.add("SymbolicRegression")
# Pkg.add("Random")

using CSV
using DataFrames
using Statistics
using Random
using Plots

import SymbolicRegression:
    Options,
    equation_search,
    calculate_pareto_frontier


# ============================================================
# 1. 基本设置
# ============================================================

const LAG_WEIGHTS = Float64[
    0.5629,
    0.1257,
    0.1178,
    0.0906,
    0.0512,
    0.0258,
    0.0147,
    0.0113,
]

@assert length(LAG_WEIGHTS) == 8
@assert isapprox(sum(LAG_WEIGHTS), 1.0; atol=1e-10)

println("Lag weights sum = ", sum(LAG_WEIGHTS))


# ============================================================
# 2. 读取和检查数据
# ============================================================

"""
    read_region(path, region)

读取单个地区数据。

要求 CSV 至少包含：
    Day
    beta
    temp
"""
function read_region(path::String, region::String)

    df = CSV.read(path, DataFrame)

    required_cols = [:Day, :beta, :temp]

    for col in required_cols
        if !(col in propertynames(df))
            error("文件 $path 缺少列: $col")
        end
    end

    # 按 Day 排序，防止原始文件顺序有问题
    sort!(df, :Day)

    # 转成 Float64
    df.beta = Float64.(df.beta)
    df.temp = Float64.(df.temp)

    # 基本检查
    if any(ismissing, df.beta)
        error("$region 的 beta 存在缺失值")
    end

    if any(ismissing, df.temp)
        error("$region 的 temp 存在缺失值")
    end

    if any(df.beta .<= 0)
        error("$region 存在 beta <= 0，无法直接进行 log(beta) 变换")
    end

    # 检查 Day 是否唯一
    if length(unique(df.Day)) != nrow(df)
        error("$region 的 Day 存在重复值")
    end

    df.region = fill(region, nrow(df))

    return df
end


# ============================================================
# 3. 构造 0-7 天滞后加权温度
# ============================================================

"""
    weighted_temperature(temp, weights)

计算：

T_eff(t) =
    w0*T(t)
  + w1*T(t-1)
  + ...
  + w7*T(t-7)

因此前 7 天无法得到完整的 T_eff，会被舍弃。
"""
function weighted_temperature(
    temp::AbstractVector{<:Real},
    weights::AbstractVector{<:Real},
)

    nlags = length(weights) - 1
    n = length(temp)

    if n <= nlags
        error("温度序列长度必须大于最大滞后阶数")
    end

    Teff = Vector{Float64}(undef, n - nlags)

    for t in (nlags + 1):n

        value = 0.0

        for k in 0:nlags
            value += weights[k + 1] * temp[t - k]
        end

        Teff[t - nlags] = value
    end

    return Teff
end


# ============================================================
# 3.5 构造14天温度变化率
# ============================================================

"""
计算14天温度变化平均速度：

DeltaT14(t) = (T(t)-T(t-14))/14

表示过去14天平均每天升温/降温幅度。
前14天无法计算。
"""
function temperature_change_14(temp::AbstractVector{<:Real})

    n = length(temp)
    lag = 30

    if n <= lag
        error("温度序列长度必须大于14天")
    end

    DeltaT14 = Vector{Float64}(undef, n-lag)

    for t in (lag+1):n
        DeltaT14[t-lag] = (temp[t] - temp[t-lag]) / lag
    end

    return DeltaT14
end


# ============================================================
# 4. 单地区预处理
# ============================================================

"""
    prepare_region(df, weights)

得到：

    Teff
    logbeta
    y = logbeta - mean(logbeta)

同时返回该地区 beta 的几何均值 baseline。
"""
function prepare_region(
    df::DataFrame,
    weights::Vector{Float64},
)

    nlags = length(weights) - 1

    temp = Float64.(df.temp)
    beta = Float64.(df.beta)

    # ---- 滞后加权温度 ----
    Teff = weighted_temperature(temp, weights)

    # ---- 14天温度变化率 ----
    DeltaT14 = temperature_change_14(temp)

    # 两个变量统一从第15天开始
    start_day = 31

    beta_valid = beta[start_day:end]
    day_valid = df.Day[start_day:end]

    Teff = Teff[(start_day - nlags):end]
    DeltaT14 = DeltaT14[(start_day - 30):end]

    @assert length(Teff) == length(beta_valid)
    @assert length(DeltaT14) == length(beta_valid)

    # ---- log beta ----
    logbeta = log.(beta_valid)
# ---- beta baseline ----
beta_mean = mean(beta_valid)

# ---- log(beta / beta_mean) ----
y = log.(beta_valid ./ beta_mean)

# 保存备用
mean_logbeta = log(beta_mean)

beta_geo = beta_mean

    # ---- 去地区 baseline ----
    y = logbeta .- mean_logbeta

    region_name = string(df.region[1])

    processed = DataFrame(
        region = fill(region_name, length(y)),
        Day = day_valid,
        temp_eff = Teff,
        temp_change14 = DeltaT14,
        beta = beta_valid,
        logbeta = logbeta,
        y = y,
    )

    return processed, mean_logbeta, beta_geo
end


# ============================================================
# 5. 读取三个地区
# ============================================================

# 修改成你自己的路径
IA_raw = read_region("C:\\Users\\XLYC\\Desktop\\流感代码修改\\AAAinfluenzacode\\SR\\IAbeta.csv", "IA")
IL_raw = read_region("C:\\Users\\XLYC\\Desktop\\流感代码修改\\AAAinfluenzacode\\SR\\ILbeta.csv", "IL")
IN_raw = read_region("C:\\Users\\XLYC\\Desktop\\流感代码修改\\AAAinfluenzacode\\SR\\INbeta.csv", "IN")


# 检查三个地区样本量
println("IA raw observations = ", nrow(IA_raw))
println("IL raw observations = ", nrow(IL_raw))
println("IN raw observations = ", nrow(IN_raw))

@assert nrow(IA_raw) == nrow(IL_raw) == nrow(IN_raw)


# ============================================================
# 6. 构造有效温度 + log baseline normalization
# ============================================================

IA, mu_IA, baseline_IA =
    prepare_region(IA_raw, LAG_WEIGHTS)

IL, mu_IL, baseline_IL =
    prepare_region(IL_raw, LAG_WEIGHTS)

IN, mu_IN, baseline_IN =
    prepare_region(IN_raw, LAG_WEIGHTS)


println()
println("Geometric mean beta:")
println("IA = ", baseline_IA)
println("IL = ", baseline_IL)
println("IN = ", baseline_IN)


# ============================================================
# 7. 合并三个地区
# ============================================================

data = vcat(IA, IL, IN)

println()
println("Total SR observations = ", nrow(data))

println()
println("Target mean by region:")
println(
    combine(
        groupby(data, :region),
        :y => mean => :mean_y,
        :y => std => :sd_y,
    ),
)


# ============================================================
# 8. 时间序列 train/test split
# ============================================================
#
# 不要随机 split。
#
# 每一个地区：
# 前 80% -> training
# 后 20% -> testing
#

function temporal_split(
    df::DataFrame;
    train_fraction::Float64 = 0.9,
)

    n = nrow(df)

    cutoff = floor(Int, train_fraction * n)

    train = df[1:cutoff, :]
    test = df[(cutoff + 1):end, :]

    return train, test
end


IA_train, IA_test = temporal_split(IA)
IL_train, IL_test = temporal_split(IL)
IN_train, IN_test = temporal_split(IN)

train_data = vcat(
    IA_train,
    IL_train,
    IN_train,
)

test_data = vcat(
    IA_test,
    IL_test,
    IN_test,
)


# ============================================================
# 9. 转成 SymbolicRegression.jl 输入格式
# ============================================================
#
# X 必须是:
#
#     [features × observations]
#
# 我们只有一个 predictor：T_eff
#

X_train = hcat(
    Float64.(train_data.temp_eff),
    Float64.(train_data.temp_change14),
)'


y_train = Float64.(train_data.y)

X_test = hcat(
    Float64.(test_data.temp_eff),
    Float64.(test_data.temp_change14),
)'


y_test = Float64.(test_data.y)


println()
println("Training observations = ", length(y_train))
println("Testing observations  = ", length(y_test))


# ============================================================
# 10. 符号回归 operator 设置
# ============================================================
#
# 第一轮建议 operator 尽可能保守。
#
# 原因：
#
# 1. 我们的目标是科学解释，不是单纯预测。
# 2. 温度只有一个自变量。
# 3. 太多 operator 很容易产生复杂但不可解释的函数。
#
# 先搜索：
#
#     +  -  *  /
#
# 即可自动产生：
#
# T
# T^2
# T^3
# 1/T
# a + b*T
# a + b*T + c*T^2
# ...
#
square(x)=x^2
cube(x)=x^3
options = Options(

    binary_operators = [
        +,
        -,
        *,
        /
    ],

    unary_operators = [
        exp,
        log,
        sqrt,
        cube,
        square   
    ],

    populations = 50,

    population_size =200,

    maxsize = 30,

    batching = false,
)


# ============================================================
# 11. 固定随机种子
# ============================================================

Random.seed!(222)


# ============================================================
# 12. 联合 Symbolic Regression
# ============================================================
#
# 注意：
#
# IA / IL / IN 已经合并。
#
# 因此这里寻找的只有一个函数：
#
#       y = g(T_eff)
#
# 不存在 g_IA / g_IL / g_IN。
#

hall_of_fame = equation_search(

    X_train,
    y_train;

    niterations = 100,

    options = options,

    variable_names = ["T_eff", "DeltaT14"],

    parallelism = :multithreading,
)


# ============================================================
# 13. 提取 Pareto frontier
# ============================================================

pareto = calculate_pareto_frontier(hall_of_fame)

println()
println("========================================")
println("Pareto-front equations")
println("========================================")

for (j, member) in enumerate(pareto)

    println()
    println("Model ", j)
    println("Expression = ", member.tree)
    println("Training loss = ", member.loss)

end


# ============================================================
# 14. 模型评价函数
# ============================================================

function metrics(
    ytrue::AbstractVector,
    ypred::AbstractVector,
)

    residual = ytrue .- ypred

    mse = mean(residual .^ 2)

    rmse = sqrt(mse)

    mae = mean(abs.(residual))

    sst = sum(
        (ytrue .- mean(ytrue)) .^ 2
    )

    ssr = sum(residual .^ 2)

    r2 = 1 - ssr / sst

    return (
        RMSE = rmse,
        MAE = mae,
        R2 = r2,
    )
end



# ============================================================
# 14.5 地区响应强度估计
# ============================================================
#
# 新模型：
#
#     y_region(t) = k_region * g(T_eff, DeltaT14)
#
# 其中：
#
#     g(.)     三个地区共享温度响应结构
#     k_region 三个地区独立估计的响应强度
#
# 对每一个候选 symbolic function，
# 先计算 g，再用最小二乘估计 k。
#

function estimate_region_k(
    gvalue::AbstractVector,
    yvalue::AbstractVector,
)

    denom = sum(gvalue .^ 2)

    if denom < 1e-12
        return 0.0
    end

    return sum(gvalue .* yvalue) / denom
end


function estimate_all_k(
    tree,
    df::DataFrame,
)

    X = hcat(
        Float64.(df.temp_eff),
        Float64.(df.temp_change14),
    )'

    g = vec(tree(X))

    ks = Dict{String, Float64}()

    for r in unique(df.region)

        idx = df.region .== r

        ks[string(r)] = estimate_region_k(
            g[idx],
            Float64.(df.y[idx]),
        )
    end

    return ks
end


function predict_beta_scaled(
    tree,
    df::DataFrame,
    beta_geo::Float64,
    k::Float64,
)

    X = hcat(
        Float64.(df.temp_eff),
        Float64.(df.temp_change14),
    )'

    g = vec(tree(X))

    return beta_geo .* exp.(k .* g)

end



# ============================================================
# 15. 在测试集评价 Pareto-front 中所有表达式
# ============================================================

println()
println("========================================")
println("Out-of-sample performance")
println("========================================")

results = DataFrame(
    model = Int[],
    expression = String[],
    train_loss = Float64[],
    test_RMSE = Float64[],
    test_MAE = Float64[],
    test_R2 = Float64[],
)

for (j, member) in enumerate(pareto)

    tree = member.tree

    # --------------------------------------------------------
    # 新结构：
    #
    # y_region = k_region * g(T)
    #
    # k 不进入 symbolic search，
    # 而是在每个地区独立估计。
    # --------------------------------------------------------

    k_train = estimate_all_k(tree, train_data)

    X_test_local = hcat(
        Float64.(test_data.temp_eff),
        Float64.(test_data.temp_change14),
    )'

    g_test = vec(tree(X_test_local))

    yhat_test = similar(g_test)

    for i in eachindex(yhat_test)

        r = string(test_data.region[i])

        yhat_test[i] = k_train[r] * g_test[i]

    end

    m = metrics(y_test, yhat_test)

    push!(
        results,
        (
            j,
            string(tree),
            Float64(member.loss),
            m.RMSE,
            m.MAE,
            m.R2,
        ),
    )
end


sort!(results, :test_RMSE)

println(results)


# ============================================================
# 16. 保存模型比较
# ============================================================

CSV.write(
    "C:\\Users\\XLYC\\Desktop\\流感代码修改\\AAAinfluenzacode\\SR\\symbolic_regression_pareto_results.csv",
    results,
)


# ============================================================
# 17. 选一个最终模型
# ============================================================
#
# 注意：
#
# 不建议机械地选择 test_RMSE 最小的表达式。
#
# 更合理：
# 在预测性能接近时优先选择结构最简单、
# 生物学上可解释的函数。
#
# 这里先演示 RMSE 最优。
#

best_index = results.model[1]

best_member = pareto[best_index]

best_tree = best_member.tree

println()
println("========================================")
println("Selected common temperature function with region-specific intensity")
println("========================================")

println("g(T_eff) = ", best_tree)


# ============================================================
# 18. 恢复到 beta 原始尺度
# ============================================================

#
# 模型：
#
# y_i(t)
# =
# log(beta_i(t)) - mean(log(beta_i))
#
# SR:
#
# y = g(T_eff, DeltaT14)
#
# 所以：
#
# beta_i(t)
# =
# beta_geo_i * exp(g(T_eff, DeltaT14))
#


function predict_beta(
    tree,
    temp_eff::Vector{Float64},
    deltaT14::Vector{Float64},
    beta_geo::Float64,
    k_region::Float64,
)

    X = hcat(
        temp_eff,
        deltaT14
    )'


    g = vec(tree(X))


    return beta_geo .* exp.(k_region .* g)

end



# ============================================================
# 三个地区预测
# ============================================================

# 独立估计 IA / IL / IN 的温度响应强度 k
k_all = estimate_all_k(best_tree, data)

println()
println("Estimated region response intensity:")
println("k_IA = ", k_all["IA"])
println("k_IL = ", k_all["IL"])
println("k_IN = ", k_all["IN"])


IA.beta_pred = predict_beta(
    best_tree,
    Float64.(IA.temp_eff),
    Float64.(IA.temp_change14),
    baseline_IA,
    k_all["IA"],
)


IL.beta_pred = predict_beta(
    best_tree,
    Float64.(IL.temp_eff),
    Float64.(IL.temp_change14),
    baseline_IL,
    k_all["IL"],
)


IN.beta_pred = predict_beta(
    best_tree,
    Float64.(IN.temp_eff),
    Float64.(IN.temp_change14),
    baseline_IN,
    k_all["IN"],
)



# ============================================================
# 19. 原始 beta 尺度评价
# ============================================================


println()
println("========================================")
println("Performance on original beta scale")
println("========================================")


println(
    "IA: ",
    metrics(
        Float64.(IA.beta),
        Float64.(IA.beta_pred),
    ),
)


println(
    "IL: ",
    metrics(
        Float64.(IL.beta),
        Float64.(IL.beta_pred),
    ),
)


println(
    "IN: ",
    metrics(
        Float64.(IN.beta),
        Float64.(IN.beta_pred),
    ),
)



# ============================================================
# 20. 保存预测数据
# ============================================================


CSV.write(
    "C:\\Users\\XLYC\\Desktop\\流感代码修改\\AAAinfluenzacode\\SR\\IA_symbolic_prediction.csv",
    IA,
)


CSV.write(
    "C:\\Users\\XLYC\\Desktop\\流感代码修改\\AAAinfluenzacode\\SR\\IL_symbolic_prediction.csv",
    IL,
)


CSV.write(
    "C:\\Users\\XLYC\\Desktop\\流感代码修改\\AAAinfluenzacode\\SR\\IN_symbolic_prediction.csv",
    IN,
)



println()
println("Finished.")



# ============================================================
# 21. 可视化 原始 beta vs 拟合 beta
# ============================================================


println()
println("Plotting beta fitting results...")



# ---------- IA ----------

p1 = plot(
    IA.Day,
    IA.beta,
    label = "Original beta",
    linewidth = 2,
    xlabel = "Day",
    ylabel = "beta",
    title = "IA: Original vs Symbolic Regression",
)


plot!(
    p1,
    IA.Day,
    IA.beta_pred,
    label = "Fitted beta",
    linewidth = 2,
)



# ---------- IL ----------


p2 = plot(
    IL.Day,
    IL.beta,
    label = "Original beta",
    linewidth = 2,
    xlabel = "Day",
    ylabel = "beta",
    title = "IL: Original vs Symbolic Regression",
)


plot!(
    p2,
    IL.Day,
    IL.beta_pred,
    label = "Fitted beta",
    linewidth = 2,
)



# ---------- IN ----------


p3 = plot(
    IN.Day,
    IN.beta,
    label = "Original beta",
    linewidth = 2,
    xlabel = "Day",
    ylabel = "beta",
    title = "IN: Original vs Symbolic Regression",
)


plot!(
    p3,
    IN.Day,
    IN.beta_pred,
    label = "Fitted beta",
    linewidth = 2,
)



display(p1)
display(p2)
display(p3)



# ============================================================
# 22. 预测散点图
# ============================================================


scatter(
    IA.beta,
    IA.beta_pred,
    xlabel="Original beta",
    ylabel="Predicted beta",
    title="IA prediction scatter",
    label="samples",
)



scatter(
    IL.beta,
    IL.beta_pred,
    xlabel="Original beta",
    ylabel="Predicted beta",
    title="IL prediction scatter",
    label="samples",
)



scatter(
    IN.beta,
    IN.beta_pred,
    xlabel="Original beta",
    ylabel="Predicted beta",
    title="IN prediction scatter",
    label="samples",
)



# ============================================================
# 23. 检查两个变量相关性
# ============================================================


println()
println("========================================")
println("Feature correlation")
println("========================================")


println(
    "Correlation T_eff = ",
    cor(data.temp_eff, data.y)
)


println(
    "Correlation DeltaT14 = ",
    cor(data.temp_change14, data.y)
)