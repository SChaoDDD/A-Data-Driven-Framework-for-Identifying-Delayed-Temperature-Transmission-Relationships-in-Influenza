using Pkg
Pkg.activate()
using SpecialFunctions
using Lux, DiffEqFlux, DifferentialEquations, Optimization, OptimizationOptimJL, Random, Plots
using DataFrames, CSV, ComponentArrays, OptimizationOptimisers
using LinearAlgebra, Distributions, Statistics
using Zygote, SciMLSensitivity

rng = Random.default_rng()
Random.seed!(14)
const SAVE_PATH = "C:/Users/XLYC/Desktop/流感代码修改/PLOS代码/JUlia/IL/ILuncertain1/"

# ============================================================
# 1. 数据准备
# ============================================================
source_data = DataFrame(CSV.File(
    "C:/Users/XLYC/Desktop/流感代码修改/PLOS代码/matlab/ILEMA.csv"))
n = 0
m = 312
data_weekly = source_data[(n+1):n+m, 1] * 2.74
display(plot(data_weekly, label="Weekly Cases", lw=2))
println("数据长度: ", length(data_weekly))

data_acc = zeros(Float32, length(data_weekly))
data_acc[1] = data_weekly[1]
for j in 2:length(data_weekly)
    data_acc[j] = sum(data_weekly[1:j])
end

trainingdata_1 = Float32.(data_acc)
trainingdata_2 = Float32.(data_weekly)
N = 12830000.0f0

# ============================================================
# 2. 网络定义（结构固定，与原代码完全一致）
#    网络结构是共享的，只有参数初始值不同
# ============================================================
ann = Lux.Chain(
    Lux.Dense(3, 16, tanh),
    Lux.Dense(16, 16, tanh),
    Lux.Dense(16, 1, softplus)
)

# ── 关键：先用原代码完全相同的方式拿到第一组参数 ──────────────
# 此时rng已经是Random.seed!(14)之后的default_rng状态，与原代码一致
p_0, st = Lux.setup(rng, ann)
# st 是网络结构状态，所有成员共用同一个 st（无状态网络）

# ============================================================
# 3. 生成各集成成员的初始参数
#    成员1 = 原代码的p_0（保证与原始结果完全一致）
#    成员2~N = 用不同MersenneTwister seed重新初始化参数
# ============================================================
const N_ENSEMBLE     = 15
const ENSEMBLE_SEEDS = [42,6,38, 77,666, 123,347, 256, 512, 999,456,88,66,100]  # 成员2~7的seed

function init_params_from_seed(seed::Int)
    # 用独立的MersenneTwister，不污染全局rng
    local_rng = Random.MersenneTwister(seed)
    p, _ = Lux.setup(local_rng, ann)   # st结构相同，只取p
    return ComponentArray(p)
end

# 收集所有成员的初始nn参数
all_init_nn_params = Vector{Any}(undef, N_ENSEMBLE)
all_init_nn_params[1] = ComponentArray(p_0)           # 成员1 = 原代码
for i in 2:N_ENSEMBLE
    all_init_nn_params[i] = init_params_from_seed(ENSEMBLE_SEEDS[i-1])
end

# ============================================================
# 4. ODE（使用全局ann和st，所有成员共享同一网络结构）
# ============================================================
function SEIR_nn(du, u, p, t)
    S, E, I, R, C = u
    t1 = sin(t * 2.0f0 * Float32(π) / 52.0f0)
    t2 = -cos(t * 2.0f0 * Float32(π) / 52.0f0) + 1.0f0
    t3 = t / 156.0f0
    nn_input = [t1, t2, t3]
    beta_raw = ann(nn_input, p, st)[1][1]
    beta     = clamp(beta_raw, 0.01f0, 5.0f0)
    S_ = max(S, 0.0f0); E_ = max(E, 0.0f0)
    I_ = max(I, 0.0f0); R_ = max(R, 0.0f0)
    du[1] = -beta * S_ * I_ / N + 0.0078f0 * R_
    du[2] =  beta * S_ * I_ / N - 3.5f0 * E_
    du[3] =  3.5f0 * E_ - I_
    du[4] =  I_ - 0.0078f0 * R_
    du[5] =  3.5f0 * E_
end

u0     = Float32[N - 9, 2, 7, 0, 7]
tspan  = (0.0f0, 311.0f0)
tsteps = range(tspan[1], tspan[2], length=312)

# prob_neuralode用p_0初始化，remake时会替换p，所以初始p无关紧要
prob_neuralode = ODEProblem(SEIR_nn, u0, tspan, ComponentArray(p_0))

const SENSE_ALG = InterpolatingAdjoint(autojacvec=ZygoteVJP())

# ============================================================
# 5. 预测函数（与原代码完全一致，所有成员共用）
# ============================================================
function predict_neuralode(θ)
    prob = remake(prob_neuralode, p=θ)
    sol  = solve(prob, Tsit5();
                 saveat=tsteps, abstol=1e-6, reltol=1e-4,
                 maxiters=100_000, verbose=false, sensealg=SENSE_ALG)
    if !SciMLBase.successful_retcode(sol) ||
       size(Array(sol), 2) != length(tsteps)
        return nothing
    end
    sol_arr = Array(sol)
    C = max.(sol_arr[5, :], 1.0f0)
    D = vcat(C[1], max.(diff(C), 1.0f0))
    return C, D
end

# 测试（与原代码一致）
println("测试初始预测...")
test_result = predict_neuralode(ComponentArray(p_0))
if test_result !== nothing
    println("✅ 初始预测成功！C[end] = $(test_result[1][end])")
else
    error("❌ 初始预测失败，检查初始条件")
end

# ============================================================
# 6. 损失函数（与原代码完全一致）
# ============================================================
function negbin_loglik(y::Real, μ::Real, r::Real)
    μ_d = max(Float64(μ), 1e-8)
    r_d = max(Float64(r), 1e-8)
    y_d = max(Float64(y), 0.0)
    ll  = loggamma(y_d + r_d) - loggamma(r_d) - loggamma(y_d + 1.0) +
          r_d * log(r_d / (r_d + μ_d)) + y_d * log(μ_d / (r_d + μ_d))
    return Float32(isfinite(ll) ? ll : -1f10)
end

function loss_negbin(p_full)
    r      = exp(Float32(p_full.log_r))
    result = predict_neuralode(p_full.nn_params)
    result === nothing && return 1.0f10
    _, pred_D = result
    nll = 0.0f0
    @inbounds for i in eachindex(trainingdata_2)
        nll -= negbin_loglik(trainingdata_2[i], pred_D[i], r)
    end
    return isfinite(nll) ? nll : 1.0f10
end

# 验证成员1的初始NLL与原代码一致
p_check = ComponentArray(nn_params=all_init_nn_params[1], log_r=Float32(log(10.0)))
println("成员#1 初始 NegBin NLL = ", loss_negbin(p_check), "  （应与原代码一致）")

# ============================================================
# 7. Beta计算函数（与原代码一致，共用全局ann/st）
# ============================================================
function compute_daily_beta(θ, days=312*7)
    betas = Float64[]
    for d in 0:(days-1)
        t  = d / 7.0
        t1 = Float32(sin(t * 2π / 52))
        t2 = Float32(-cos(t * 2π / 52) + 1)
        t3 = Float32(t / 156)
        β  = clamp(Float64(ann([t1, t2, t3], θ, st)[1][1]), 0.01, 5.0)
        push!(betas, β)
    end
    return betas
end

# ============================================================
# 8. 单成员训练函数
# ============================================================
function train_single_member(id::Int, init_nn_params)
    println("\n" * "="^60)
    println("🔄 训练集成成员 #$id")
    println("="^60)

    p_init = ComponentArray(
        nn_params = init_nn_params,
        log_r     = Float32(log(10.0))
    )
    println("  初始 NLL = $(loss_negbin(p_init))")

    iter_cnt = Ref(0)
    cb = function(state, l)
        isfinite(l) && (iter_cnt[] += 1)
        iter_cnt[] % 500 == 0 &&
            println("  [#$id] Iter $(iter_cnt[]): NLL = $(round(l, digits=2))")
        return false
    end

    adtype  = Optimization.AutoZygote()
    optf    = Optimization.OptimizationFunction((x, _) -> loss_negbin(x), adtype)
    optprob = Optimization.OptimizationProblem(optf, p_init)

    # —— ADAM ——
    res_adam = Optimization.solve(optprob,
        OptimizationOptimisers.ADAM(0.001), callback=cb, maxiters=3500)
    nll_after_adam = loss_negbin(res_adam.u)
    println("  [#$id] ADAM  NLL = $(round(nll_after_adam, digits=2))")

    # —— L-BFGS（以ADAM结果为起点，allow_f_increases=false保证不变差）——
    pfinal = try
        optprob2 = remake(optprob, u0=res_adam.u)
        res_lb   = Optimization.solve(optprob2, Optim.LBFGS();
                       callback          = cb,
                       allow_f_increases = false,
                       maxiters          = 3000)
        nll_lb = loss_negbin(res_lb.u)
        println("  [#$id] LBFGS NLL = $(round(nll_lb, digits=2))")
        # 如果LBFGS结果比ADAM更差，放弃LBFGS结果
        nll_lb <= nll_after_adam ? res_lb.u : res_adam.u
    catch e
        println("  [#$id] LBFGS 异常 ($e) → 保留 ADAM 结果")
        res_adam.u
    end

    final_nll = loss_negbin(pfinal)
    r_val     = exp(Float64(pfinal.log_r))
    println("  [#$id] ✅ 完成  NLL = $(round(final_nll, digits=2))  " *
            "r = $(round(r_val, digits=3))")

    return (id       = id,
            params   = pfinal,
            nll      = final_nll,
            r        = r_val)
end

# ============================================================
# 9. 深度集成：依次训练所有成员
# ============================================================
println("\n🚀 开始深度集成训练（$N_ENSEMBLE 个成员）")

members = [train_single_member(i, all_init_nn_params[i]) for i in 1:N_ENSEMBLE]

# ============================================================
# 10. 筛选：最优基准 + 剔除差模型
# ============================================================
nll_values = [m.nll for m in members]

println("\n" * "="^60)
println("📋 各成员 NLL 汇总")
println("-"^40)
for m in members
    println("  成员#$(m.id)  NLL = $(round(m.nll, digits=2))  " *
            "r = $(round(m.r, digits=3))")
end

# —— 最优基准 ——
best_idx = argmin(nll_values)
best_nll = nll_values[best_idx]
println("\n⭐ 最佳基准: 成员#$(members[best_idx].id)  " *
        "NLL = $(round(best_nll, digits=2))")

# —— 剔除阈值：median + 1.5 × IQR（箱线图法）——
finite_nlls = filter(x -> x < 1.0f9, nll_values)
nll_cutoff  = if length(finite_nlls) >= 3
    q1, q3  = quantile(finite_nlls, 0.25), quantile(finite_nlls, 0.75)
    iqr      = q3 - q1
    # IQR过小时用best_nll的2%作为最小间隔
    median(finite_nlls) + 1.5 * max(iqr, 0.02 * abs(best_nll))
else
    best_nll + abs(best_nll) * 0.10
end
println("  剔除阈值 NLL > $(round(nll_cutoff, digits=2))")

keep_mask = [(m.nll < 1.0f9) && (m.nll <= nll_cutoff)&&(m.r>3) for m in members]

println("\n筛选结果:")
for (i, m) in enumerate(members)
    tag  = keep_mask[i] ? "✅ 保留" : "❌ 剔除"
    star = (i == best_idx) ? " ⭐基准" : ""
    println("  成员#$(m.id)  NLL=$(round(m.nll, digits=2))  $tag$star")
end

kept_indices = findall(keep_mask)
n_kept       = length(kept_indices)

if n_kept < 2
    println("⚠️ 保留成员不足2个，放宽为所有有限NLL成员")
    keep_mask    = [m.nll < 1.0f9 for m in members]
    kept_indices = findall(keep_mask)
    n_kept       = length(kept_indices)
end
println("\n最终保留 $n_kept / $N_ENSEMBLE 个成员")

# ============================================================
# 11. 集成预测 + 不确定性
# ============================================================
println("\n计算集成预测和不确定性...")
n_weeks = 312
n_days  = 312 * 7

# 新增病例矩阵
all_pred_D = fill(NaN, n_weeks, n_kept)
all_pred_C = fill(NaN, n_weeks, n_kept)

for (col, idx) in enumerate(kept_indices)
    res = predict_neuralode(members[idx].params.nn_params)
    if res !== nothing
        C_i, D_i = res
        all_pred_C[:, col] = Float64.(C_i)
        all_pred_D[:, col] = Float64.(D_i)
    else
        println("⚠️ 成员#$(members[idx].id) 集成预测失败")
    end
end

# 去除NaN列
valid_mask   = [!any(isnan.(all_pred_D[:, c])) for c in 1:n_kept]
all_pred_D   = all_pred_D[:, valid_mask]
all_pred_C   = all_pred_C[:, valid_mask]
final_kept   = kept_indices[valid_mask]
n_final      = length(final_kept)
println("有效Ensemble成员数: $n_final")

# 统计量
ens_D_mean  = vec(mean(all_pred_D, dims=2))
ens_D_std   = vec(std(all_pred_D,  dims=2))
ens_D_lower = ens_D_mean .- 1.96 .* ens_D_std
ens_D_upper = ens_D_mean .+ 1.96 .* ens_D_std
ens_D_q025  = vec(mapslices(x -> quantile(x, 0.025), all_pred_D, dims=2))
ens_D_q975  = vec(mapslices(x -> quantile(x, 0.975), all_pred_D, dims=2))

# Beta矩阵
all_betas = fill(NaN, n_days, n_final)
for (col, idx) in enumerate(final_kept)
    all_betas[:, col] = compute_daily_beta(members[idx].params.nn_params, n_days)
end

ens_beta_mean  = vec(mean(all_betas, dims=2))
ens_beta_std   = vec(std(all_betas,  dims=2))
ens_beta_lower = max.(ens_beta_mean .- 1.96 .* ens_beta_std, 0.01)
ens_beta_upper = min.(ens_beta_mean .+ 1.96 .* ens_beta_std, 5.0)
ens_beta_q025  = max.(vec(mapslices(x -> quantile(x, 0.025), all_betas, dims=2)), 0.01)
ens_beta_q975  = min.(vec(mapslices(x -> quantile(x, 0.975), all_betas, dims=2)), 5.0)

# 最佳基准模型的预测
best_m          = members[best_idx]
best_C, best_D  = predict_neuralode(best_m.params.nn_params)
best_beta       = compute_daily_beta(best_m.params.nn_params, n_days)

# ============================================================
# 12. 保存 CSV
# ============================================================
println("\n保存 CSV...")

df_beta = DataFrame(
    Day           = 1:n_days,
    Beta_Best     = best_beta,
    Beta_EnsMean  = ens_beta_mean,
    Beta_EnsStd   = ens_beta_std,
    Beta_Lower95  = ens_beta_lower,
    Beta_Upper95  = ens_beta_upper,
    Beta_Q025     = ens_beta_q025,
    Beta_Q975     = ens_beta_q975
)
CSV.write(SAVE_PATH * "daily_beta_ensemble.csv", df_beta)
println("✅ daily_beta_ensemble.csv  ($(nrow(df_beta)) 行)")

df_weekly = DataFrame(
    Week          = 1:n_weeks,
    Observed      = Float64.(trainingdata_2),
    Pred_Best     = Float64.(best_D),
    Pred_EnsMean  = ens_D_mean,
    Pred_EnsStd   = ens_D_std,
    Pred_Lower95  = ens_D_lower,
    Pred_Upper95  = ens_D_upper,
    Pred_Q025     = ens_D_q025,
    Pred_Q975     = ens_D_q975
)
CSV.write(SAVE_PATH * "weekly_cases_ensemble.csv", df_weekly)
println("✅ weekly_cases_ensemble.csv  ($(nrow(df_weekly)) 行)")

# ============================================================
# 13. 可视化
# ============================================================

# —— Figure 1: 每周新增病例 ——
plt1 = plot(1:n_weeks, Float64.(trainingdata_2);
    label="Observed Data", lw=2, color=:blue,
    marker=:circle, markersize=2, alpha=0.8,
    xlabel="Week", ylabel="New Cases",
    title="Weekly New Cases — Best Model + Ensemble CI  ($n_final/$N_ENSEMBLE members)",
    size=(1100, 500), legend=:topright)

plot!(plt1, 1:n_weeks, ens_D_q025;
    fillrange=ens_D_q975, fillalpha=0.22, fillcolor=:orange,
    label="95% Ensemble CI", linealpha=0)

for col in 1:n_final
    plot!(plt1, 1:n_weeks, all_pred_D[:, col];
        label=(col == 1 ? "Ensemble Members" : ""),
        lw=0.8, color=:gray, alpha=0.45)
end

plot!(plt1, 1:n_weeks, ens_D_mean;
    label="Ensemble Mean", lw=2, color=:darkorange, linestyle=:dot)
plot!(plt1, 1:n_weeks, Float64.(best_D);
    label="Best Model (#$(best_m.id))", lw=2.5, color=:red, linestyle=:dash)

display(plt1)
savefig(plt1, SAVE_PATH * "weekly_cases_ensemble_CI.png")

# —— Figure 2: 每日 Beta ——
plt2 = plot(1:n_days, best_beta;
    label="Best Model β (#$(best_m.id))", lw=2, color=:green,
    xlabel="Day", ylabel="β",
    title="β(t) — Best Model + Ensemble CI  ($n_final/$N_ENSEMBLE members)",
    size=(1200, 450), legend=:topright)

plot!(plt2, 1:n_days, ens_beta_q025;
    fillrange=ens_beta_q975, fillalpha=0.22, fillcolor=:green,
    label="95% Ensemble CI", linealpha=0)

for col in 1:n_final
    plot!(plt2, 1:n_days, all_betas[:, col];
        label=(col == 1 ? "Ensemble Members" : ""),
        lw=0.5, color=:gray, alpha=0.35)
end

plot!(plt2, 1:n_days, ens_beta_mean;
    label="Ensemble Mean", lw=1.5, color=:darkgreen, linestyle=:dot)

for yr in 1:6
    vline!(plt2, [yr*365]; label="", color=:gray, linestyle=:dash, alpha=0.5)
end

display(plt2)
savefig(plt2, SAVE_PATH * "daily_beta_ensemble_CI.png")

# —— Figure 3: 累积病例 ——
plt3 = plot(1:n_weeks, Float64.(trainingdata_1);
    label="Observed Cumulative", lw=2, color=:blue,
    xlabel="Week", ylabel="Cumulative Cases",
    title="Cumulative Cases — Best Model (#$(best_m.id))", size=(1000, 400))
plot!(plt3, 1:n_weeks, Float64.(best_C);
    label="Best Model", lw=2, color=:red, linestyle=:dash)
display(plt3)
savefig(plt3, SAVE_PATH * "cumulative_ensemble.png")

# —— Figure 4: 残差分析 ——
residuals = Float64.(trainingdata_2) .- Float64.(best_D)
plt4 = plot(layout=(1, 2), size=(1200, 420))
plot!(plt4[1], 1:n_weeks, residuals;
    label="Residuals", lw=1.2, color=:darkblue,
    xlabel="Week", ylabel="Residuals",
    title="Residuals — Best Model (#$(best_m.id))")
hline!(plt4[1], [0]; label="", color=:gray, linestyle=:dash)
histogram!(plt4[2], residuals;
    bins=30, label="Distribution", color=:steelblue, alpha=0.7,
    xlabel="Residuals", ylabel="Frequency", title="Residual Histogram")
display(plt4)
savefig(plt4, SAVE_PATH * "residuals_ensemble.png")

# ============================================================
# 14. 拟合指标
# ============================================================
obs = Float64.(trainingdata_2)

println("\n" * "="^60)
println("📊 Goodness of Fit Metrics")
println("="^60)

println("\n【最佳基准模型 #$(best_m.id)】")
println("  NLL  = $(round(best_m.nll, digits=2))")
println("  r    = $(round(best_m.r,   digits=3))")
mse_b = mean((obs .- Float64.(best_D)).^2)
mae_b = mean(abs.(obs .- Float64.(best_D)))
mape_b = mean(abs.(obs .- Float64.(best_D)) ./ max.(obs, 1.0)) * 100
println("  MSE  = $(round(mse_b,       digits=2))")
println("  RMSE = $(round(sqrt(mse_b), digits=2))")
println("  MAE  = $(round(mae_b,       digits=2))")
println("  MAPE = $(round(mape_b,      digits=2))%")

println("\n【Ensemble 不确定性（$n_final 个保留成员）】")
in_q = sum((obs .>= ens_D_q025) .& (obs .<= ens_D_q975))
in_s = sum((obs .>= ens_D_lower) .& (obs .<= ens_D_upper))
println("  95% CI Coverage (quantile): $(round(in_q/n_weeks*100, digits=1))%")
println("  95% CI Coverage (±1.96σ):   $(round(in_s/n_weeks*100, digits=1))%")
println("  Avg CI Width:               $(round(mean(ens_D_q975 .- ens_D_q025), digits=2))")

println("\n【Ensemble 均值】")
mse_e  = mean((obs .- ens_D_mean).^2)
mae_e  = mean(abs.(obs .- ens_D_mean))
mape_e = mean(abs.(obs .- ens_D_mean) ./ max.(obs, 1.0)) * 100
println("  MSE  = $(round(mse_e,       digits=2))")
println("  RMSE = $(round(sqrt(mse_e), digits=2))")
println("  MAE  = $(round(mae_e,       digits=2))")
println("  MAPE = $(round(mape_e,      digits=2))%")
println("="^60)