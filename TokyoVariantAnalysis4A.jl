module TokyoVariantAnalysis4A
using Logging
using CSV
using DataFrames
using Dates
using LinearAlgebra: dot
using SpecialFunctions: beta_inc
using LsqFit: curve_fit, coef

using ..TokyoVariantAnalysis4

export likelihoods_df, likelihood_max

@info "module", @__MODULE__

@info "========"
@info "特定変異オッズ 最尤値探索"

# モデル
lh_model(t; α, β) = exp(α + β * t)

# 累積分布関数
# 不完全ベータ関数により二項分布の累積分布を連続値関数として表す
function Φ(z; n, zm)
    @assert z ≥ 0.0
    @assert n > 0
    @assert zm ≥ 0.0
    q = zm / (zm + 1.0)
    y = n * z / (z + 1.0)
    return beta_inc(n - y + 0.5, y + 0.5, 1.0 - q)[1]
end
# 確率密度関数（数値微分による）
const h_float64 = 2^-12
function φ(z; n, zm)
    @assert z ≥ 0.0
    @assert n > 0
    @assert zm ≥ 0.0
    h::Float64 = max(z, 1.0) * h_float64
    if (h ≤ z)
        return (Φ(z + h; n=n, zm=zm) - Φ(z - h; n=n, zm=zm)) / (2.0 * h)
    else
        return (Φ(z + h; n=n, zm=zm) - Φ(z; n=n, zm=zm)) / h
    end
end

const lh_ps_niter       =    5
const lh_ps_hwidth_α    =  5.0
const lh_ps_hwidth_β    =  1.0
const lh_ps_nslice_α    =  200
const lh_ps_nslice_β    =  200
const lh_nslice_α       =  500
const lh_nslice_β       = 1000
const lh_logL_threshold = -6.0
@info "予備探索回数", lh_ps_niter
@info "予備探索パラメータ初期探索幅", lh_ps_hwidth_α, lh_ps_hwidth_β
@info "予備探索パラメータ・スライス数", lh_ps_nslice_α, lh_ps_nslice_β
@info "本探索パラメータ・スライス数", lh_nslice_α, lh_nslice_β
@info "打ち切り尤度閾値", lh_logL_threshold

function search_odd_likelihoods(df, varname, period)::DataFrame
    oddname = varname*"_odds"
    (sdate_start, sdate_end) = period
    @info oddname
    @info sdate_start, sdate_end
    @info "--------"
    @info "データ点収集"
    @info "対象行のみを取り出した DataFrame を作成する"
    sample_df = filter(
        (row) -> (begin
            var_n = row[varname]
            odds = row[oddname]
            (sdate_start ≤ row.date_end) && (row.date_start ≤ sdate_end) &&
            !ismissing(var_n) && !ismissing(odds) && !isinf(odds) && !isnan(odds)
        end),
        df
    )
    @info size(sample_df)
    @info "--------"
    @info "サンプルサイズ"
    sample_size = nrow(sample_df)
    @info sample_size
    @info "--------"
    @info "基準日算出"
    date_epoch = round(dot(sample_df.date_value, sample_df.n) / sum(sample_df.n))
    @info date_epoch
    (t_start, t_end) = date_to_value.(period) .- date_epoch
    @info t_start, t_end
    sample_df.t = sample_df.date_value .- date_epoch
    @info "--------"
    @info "非線形フィッティングの予備的実行"
    @info "初期パラメーターを定めるため非線形最小二乗フィッティング"
    lsq_α = nothing
    lsq_β = nothing
    for once in 1:1
        cils = sample_df[!,oddname*"_cil"]
        cius = sample_df[!,oddname*"_ciu"]
        ws = 1.0 ./ (cius .- cils).^2
        βs0 = [0.0, 0.0]
        cf = curve_fit(
            (t, βs) -> lh_model.(t; α=βs[1], β=βs[2]),
            sample_df.t, sample_df[!,oddname], ws, βs0
        )
        (lsq_α, lsq_β) = coef(cf)
    end
    @info lsq_α, lsq_β
    @info "--------"
    @info "尤度予備探索"
    @info "計算すべきパラメーターの範囲を絞り込む"
    mn_α_range = nothing
    mn_β_range = nothing
    for once in 1:1
        α_range = range(lsq_α - lh_ps_hwidth_α, lsq_α + lh_ps_hwidth_α; length=lh_ps_nslice_α)
        β_range = range(lsq_β - lh_ps_hwidth_β, lsq_β + lh_ps_hwidth_β; length=lh_ps_nslice_β)
        for iter in 1:lh_ps_niter
            @info "予備探索 $(iter) 回目"
            @info α_range
            @info β_range
            ps_df = DataFrame(α=Float64[], β=Float64[], logL=Float64[])
            logL_max::Union{Float64,Nothing} = nothing
            for α in α_range
                for β in β_range
                    logL = 0.0
                    for row in eachrow(sample_df)
                        zm = lh_model(row.t; α=α, β=β)
                        logL += log(φ(row[oddname]; n=row.n, zm=zm))
                    end
                    if isnothing(logL_max) || (logL_max < logL)
                        logL_max = logL
                    end
                    if logL_max + lh_logL_threshold ≤ logL
                        push!(ps_df, [α, β, logL])
                    end
                end
            end
            @info logL_max
            α_min = last(α_range)
            α_max = first(α_range)
            β_min = last(β_range)
            β_max = first(β_range)
            for row in eachrow(ps_df)
                if logL_max + lh_logL_threshold ≤ row.logL
                    if row.α < α_min; α_min = row.α end
                    if α_max < row.α; α_max = row.α end
                    if row.β < β_min; β_min = row.β end
                    if β_max < row.β; β_max = row.β end
                end
            end
            α_torr = 1.5 * step(α_range)
            β_torr = 1.5 * step(β_range)
            α_range = range(α_min - α_torr, α_max + α_torr; length=lh_ps_nslice_α)
            β_range = range(β_min - β_torr, β_max + β_torr; length=lh_ps_nslice_β)
        end
        mn_α_range = range(first(α_range), last(α_range); length=lh_nslice_α)
        mn_β_range = range(first(β_range), last(β_range); length=lh_nslice_β)
    end
    @info "--------"
    @info "尤度本探索"
    @info "パラメーターの範囲を分割し、各パラメーターにおける尤度を計算"
    @info mn_α_range
    @info mn_β_range
    lh_df::Union{DataFrame,Nothing} = nothing
    for once in 1:1
        mn_df = DataFrame(α=Float64[], β=Float64[], logL=Float64[])
        logL_max::Union{Float64,Nothing} = nothing
        for (iα, α) in enumerate(mn_α_range)
            if iα % 10 == 0
                print(stderr, "\r$(round(Int, 100.0 * (iα / lh_nslice_α)))%")
            end 
            for β in mn_β_range
                logL = 0.0
                for row in eachrow(sample_df)
                    zm = lh_model(row.t; α=α, β=β)
                    logL += log(φ(row[oddname]; n=row.n, zm=zm))
                end
                if isnothing(logL_max) || (logL_max < logL)
                    logL_max = logL
                end
                if logL_max + lh_logL_threshold ≤ logL
                    push!(mn_df, [α, β, logL])
                end
            end
        end
        print(stderr, "\r100%\n")
        lh_df = filter((row) -> logL_max + lh_logL_threshold ≤ row.logL, mn_df)
        for i in 1:nrow(lh_df)
            lh_df.α[i] = lh_df.α[i] - lh_df.β[i] * date_epoch
        end
        sort!(lh_df, :logL; rev=true)
        @info size(lh_df)
    end
    return lh_df
end

lh_varname = "L452R"
lh_sdate = Date("2021-06-01")
lh_edate = Date("2021-08-31")
lh_sdate_value = date_to_value(lh_sdate)
lh_edate_value = date_to_value(lh_edate)

likelihoods_df = search_odd_likelihoods(variants_df, lh_varname, (lh_sdate, lh_edate))
@info "size(likelihoods_df)", size(likelihoods_df)
likelihood_max = likelihoods_df[1,:]
@info "likelihood_max", likelihood_max

@info "--------"
@info "尤度データ書き出し"
const lh_output_csv_filename = "tokyo_variant_analysis4_likelihoods.csv"
@info "lh_output_csv_filename", lh_output_csv_filename
CSV.write(lh_output_csv_filename, likelihoods_df; missingstring="NA")

end #module
