module TokyoConfirmed5
using Logging
using CSV
using DataFrames
using Dates
using Statistics: mean, stdm
using LsqFit: curve_fit
#using SingularSpectrumAnalysis: analyze
using RCall
using Printf
using Plots
using Plots.PlotMeasures

export date_to_value, value_to_date

@info "module", @__MODULE__

function has_name(df::DataFrame, name::String)
    return !isnothing(findfirst(isequal(name), names(df)))
end

@info "--------"
@info "作業ディレクトリ"
@info pwd()

@info "--------"
@info "基準日・開始日"
const epoch_date = Date("2021-01-01")
@info "epoch_date", epoch_date
#
date_to_value(d::Date) = Dates.value(d - epoch_date)
value_to_date(v)::Date = epoch_date + Day(v)
#
const start_date = Date("2021-01-01")
@info "start_date", start_date

@info "========"
@info "データ読み込み・変換"

function load_df(csv_filepath)::DataFrame
    @info csv_filepath

    @info "--------"
    @info "CSV 読込み"
    csv_df = CSV.read(csv_filepath, DataFrame; missingstring="NA")
    @info size(csv_df)

    @info "--------"
    @info "データフレーム作成"
    df = DataFrame()
    @info "--------"
    @info "日付 (date_start, date_end, days, date_error, date_value, date_mid)"
    @assert has_name(csv_df, "date")
    dates = csv_df.date
    df.date_start  = dates
    df.date_end    = dates
    df.days        = date_to_value.(df.date_end) .+ 1.0 .- date_to_value.(df.date_start)
    df.date_error  = df.days ./ 2.0
    df.date_value  = date_to_value.(df.date_start) .+ df.date_error
    df.date_mid    = value_to_date.(round.(df.date_value .- 0.5))
    @assert all(df.date_start .≤ df.date_mid .≤ df.date_end)
    @info "df.date_start[end]", df.date_start[end]
    @info "df.date_mid[end]", df.date_mid[end]
    @info "df.date_end[end]", df.date_end[end]
    @info size(df)
    n_rows = nrow(df)
    @info "--------"
    @info "曜日・祝日・休日 (holiday, holiday_next)"
    df.day_of_week = dayofweek.(df.date_mid)
    has_name(csv_df, "h")
    holidays = Vector{Bool}(undef, nrow(df))
    for i in 1:n_rows
        w = df.day_of_week[i]
        h = csv_df.h[i]
        holidays[i] = ((h == "祝") || (h == "休")) && (w != 6) && (w != 7)
    end
    df.holiday = holidays
    df.holiday_next = ([false; holidays])[1:end-1]
    @info size(df)
    @info "--------"
    @info "感染確認者数のコピー (count)"
    @assert has_name(csv_df, "count")
    df.count = csv_df.count
    @info size(df)
    @info "--------"
    @info "開始日以降への制限"
    df2 = filter((row) -> (row.date_start ≥ start_date), df)
    @info size(df2)
    return df2
end

const tokyo_confirmed_csv_filename = "東京感染確認者4.csv"
confirmed_df = load_df(tokyo_confirmed_csv_filename)

@info "========"
@info "時系列解析"

const ci95s = 1.95996398454005

const ssa_window = 7
@info "ssa_window", ssa_window

function ssanalysis!(df::DataFrame)
    n_rows = nrow(df)
    @info "--------"
    @info "予備特異スペクトル解析"
    @info "カウントの対数をトレンド, 周期成分, 残差への分解"
    df.lcount = log.(df.count .+ 1) # -Inf を防ぐため 1 だけ加える
    psa_x = df.lcount
    psa_stl = rcopy(R"""
    psa.ts <- ts($psa_x, frequency=$ssa_window, start=1)
    psa.stl <- stl(psa.ts, s.window="periodic")
    """)
    psa_rts = psa_stl[:time_series]
    df.psa_lseasonal = psa_rts[:,1]
    df.psa_ltrend    = psa_rts[:,2]
    df.psa_lresidual = psa_rts[:,3]
    df.psa_lpeculiar = df.psa_lseasonal .+ df.psa_lresidual
    @info size(df)
    @info "曜日別周期成分値"
    psa_lseasonal_dows = zeros(Float64, 7)
    for i in n_rows-6:n_rows
        w = df.day_of_week[i]
        psa_lseasonal_dows[w] = df.psa_lseasonal[i]
    end
    @info psa_lseasonal_dows
    @info "--------"
    @info "予備休日補正 (pha_lcount)"
    @info "休日翌日の補正量を推定"
    pha_df = filter((row) -> row.holiday_next, df)  # 休日翌日行抽出
    pha_n_rows = nrow(pha_df)
    @info size(pha_df)
    pha_ldecs = zeros(Float64, pha_n_rows)
    for i in 1:pha_n_rows
        w = pha_df.day_of_week[i]
        e = pha_df.psa_ltrend[i] + psa_lseasonal_dows[w]
        pha_ldecs[i] = pha_df.lcount[i] - e
    end
    pha_ldecm = mean(pha_ldecs)
    @info "pha_ldecm", pha_ldecm, exp(pha_ldecm)
    @info "休日翌日の値を補正"
    pha_lcounts = Vector{Float64}(df.lcount)
    for i in 1:n_rows
        if df.holiday_next[i]
            pha_lcounts[i] -= pha_ldecm
        end
    end
    df.pha_lcount = pha_lcounts
    @info size(df)
    @info "--------"
    @info "特異スペクトル解析 (ltrend, lseasonal)"
    sa_x = df.pha_lcount
    sa_stl = rcopy(R"""
    sa.ts <- ts($sa_x, frequency=$ssa_window, start=1)
    sa.stl <- stl(sa.ts, s.window=$ssa_window, l.degree=1, t.degree=1, s.degree=1, robust=T)
    plot(sa.stl)
    sa.stl
    """)
    sa_rts = sa_stl[:time_series]
    @info size(sa_rts)
    df.lseasonal  = sa_rts[:,1]
    df.ltrend     = sa_rts[:,2]
    df.lresidual  = sa_rts[:,3]
    df.lpeculiear = df.lseasonal .+ df.lresidual
    @info size(df)
    @info "--------"
    @info "曜日別周期成分値"
    lsdow_df = DataFrame(day_of_week = 1:7)
    lsdows_vs = [Float64[] for w in 1:7]
    for i in 1:n_rows
        w = df.day_of_week[i]
        push!(lsdows_vs[w], df.lseasonal[i])
    end
    lsdow_df.mean = [mean(lsdows_vs[w]) for w in 1:7]
    lsdow_df.sd   = [stdm(lsdows_vs[w], lsdow_df.mean[w]) for w in 1:7]
    lsdow_df.cil  = lsdow_df.mean .- (ci95s .* lsdow_df.sd)
    lsdow_df.ciu  = lsdow_df.mean .+ (ci95s .* lsdow_df.sd)
    lsdow_df.magn = exp.(lsdow_df.mean)
    @info lsdow_df
    @info "--------"
    @info "予備休日補正 (ha_lcount)"
    @info "休日翌日の補正量を推定"
    ha_df = filter((row) -> row.holiday_next, df)  # 休日翌日行抽出
    ha_n_rows = nrow(ha_df)
    @info size(ha_df)
    ha_ldecs = zeros(Float64, ha_n_rows)
    for i in 1:ha_n_rows
        w = ha_df.day_of_week[i]
        e = ha_df.ltrend[i] + lsdow_df.mean[w]
        ha_ldecs[i] = ha_df.lcount[i] - e
    end
    ha_ldec_mean = mean(ha_ldecs)
    ha_ldec_sd   = stdm(ha_ldecs, ha_ldec_mean)
    @info "ha_ldec_mean", ha_ldec_mean, exp(ha_ldec_mean)
    @info "ha_ldec_sd", ha_ldec_sd
    @info "休日翌日の値を補正"
    ha_lcounts = Vector{Float64}(df.lcount)
    for i in 1:n_rows
        if df.holiday_next[i]
            ha_lcounts[i] -= ha_ldec_mean
        end
    end
    df.ha_lcount = ha_lcounts
    df.ha_count  = exp.(ha_lcounts) .- 1.0
    @info size(df)
    @info "--------"
    @info "残差誤差 (ltrend_cil, ltrend_ciu)"
    lrs = filter((row) -> (!row.holiday_next), df).lresidual
    lrm = mean(lrs)
    lre = ci95s * stdm(lrs, lrm)
    @info lrm, lre
    df.ltrend_cil = df.ltrend .- (lre + lrm)
    df.ltrend_ciu = df.ltrend .+ (lre - lrm)
    @info size(df)
    @info "--------"
    @info "補正トレンド (trend, trend_cil, trend_ciu)"
    df.trend     = exp.(df.ltrend)     .- 1.0
    df.trend_cil = exp.(df.ltrend_cil) .- 1.0
    df.trend_ciu = exp.(df.ltrend_ciu) .- 1.0
    for i in 1:n_rows
        @assert df.trend[i] ≥ 0.0
    end
    @info size(df)
    @info "--------"
    @info "拡大係数 (lgrowth, lgrowth_cil, lgrowth_ciu)"
    xs  = df.trend
    ξs  = log.(xs)
    ξds = zeros(Float64, n_rows)
    for i in 1:n_rows
        if     i == 1
            ξd = ξs[i+1] - ξs[i]
        elseif i == n_rows
            ξd = ξs[i] - ξs[i-1]
        else
            ξd = (ξs[i+1] - ξs[i-1]) / 2.0
        end
        ξds[i] = ξd
    end
    df.lgrowth = ξds
    @info size(df)
    return (lsdow_df, ha_ldec_mean, ha_ldec_sd)
end

(lsdow_df, ha_ldec_mean, ha_ldec_sd) = ssanalysis!(confirmed_df)

@info "========"
@info "株別推定"
@info "--------"
@info "尤度データ取り込み"
const likelihoods_csv_filepath = "tokyo_variant_analysis4_likelihoods.csv"
@info likelihoods_csv_filepath
likelihoods_df = CSV.read(likelihoods_csv_filepath, DataFrame; missingstring="NA")
@info size(likelihoods_df)
@info "--------"
@info "最大尤度"
likelihood_max = likelihoods_df[1,:]
likelihood_max_α = likelihood_max.α
likelihood_max_β = likelihood_max.β
@info likelihood_max
#
δ_logit(t; α, β) = α + β * t
δ_odds(t; α, β) = exp(δ_logit(t; α=α, β=β))
δ_prop(t; α, β) = 1.0 / (1.0 + 1.0 / δ_odds(t; α=α, β=β))
δ_prop_maxL(t) = δ_prop(t; α=likelihood_max_α, β=likelihood_max_β)
# δ_prop_maxL(t) は、時刻 t における最尤線が示す δ 株の割合 
#
const likelihood_lags = [0, 7, 14]
function variants_estimate!(df::DataFrame)
    n_rows = nrow(df)
    @info "--------"
    @info "δ 株と非 δ 株のトレンド (trend_δ_lagXX, trend_nonδ_lagXX)"
    ts = df.date_value
    xs = df.trend
    for lag in likelihood_lags
        ps = δ_prop_maxL.(ts .+ lag)    # 各時刻での δ 株の割合（遅れ lag を見込む）
        xs_δ    = ps .* xs              # δ 株の感染確認者トレンド数
        xs_nonδ = xs .- xs_δ            # δ 株以外の感染確認者トレンド数
        lagname = @sprintf("lag%02d", lag)        
        df[!,"trend_δ_"*lagname]    = xs_δ
        df[!,"trend_nonδ_"*lagname] = xs_nonδ
    end
    @info "--------"
    @info "δ 株拡大係数 (lgrowth_δ_lagXX)"
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_δ = df[!,"trend_δ_"*lagname]
        ξs_δ  = log.(xs_δ)
        ξds_δ = zeros(Float64, n_rows)
        for i in 1:n_rows
            if     i == 1
                ξds_δ[i] = ξs_δ[i+1] - ξs_δ[i]
            elseif i == n_rows
                ξds_δ[i] = ξs_δ[i] - ξs_δ[i-1]
            else
                ξds_δ[i] = (ξs_δ[i+1] - ξs_δ[i-1]) / 2.0
            end
        end
        df[!,"lgrowth_δ_"*lagname] = ξds_δ
    end
    @info "--------"
    @info "非 δ 株拡大係数 (lgrowth_nonδ_lagXX)"
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_nonδ = df[!,"trend_nonδ_"*lagname]
        ξs_nonδ  = log.(xs_nonδ)
        ξds_nonδ = zeros(Float64, n_rows)
        for i in 1:n_rows
            if     i == 1
                ξds_nonδ[i] = ξs_nonδ[i+1] - ξs_nonδ[i]
            elseif i == n_rows
                ξds_nonδ[i] = ξs_nonδ[i] - ξs_nonδ[i-1]
            else
                ξds_nonδ[i] = (ξs_nonδ[i+1] - ξs_nonδ[i-1]) / 2.0
            end
        end
        df[!,"lgrowth_nonδ_"*lagname] = ξds_nonδ
    end
end

variants_estimate!(confirmed_df)

@info "--------"
@info "解析データ書き出し"
const output_csv_filename = "tokyo_confirmed_analysis4.csv"
@info output_csv_filename
CSV.write(output_csv_filename, confirmed_df; missingstring="NA")
@info size(confirmed_df)

@info "========"
@info "予測"

const holiday_list = [
    Date("2021-01-01"), Date("2021-01-02"), Date("2021-01-03"),
    Date("2021-01-11"),
    Date("2021-02-11"), Date("2021-02-23"),
    Date("2021-03-20"),
    Date("2021-04-29"),
    Date("2021-05-03"), Date("2021-05-04"), Date("2021-05-05"),
    Date("2021-07-22"), Date("2021-07-23"),
    Date("2021-08-08"), Date("2021-08-09"),
    Date("2021-08-13"), Date("2021-08-14"), Date("2021-08-15"), Date("2021-08-16"),
    Date("2021-09-20"), Date("2021-09-23"),
    Date("2021-11-03"), Date("2021-11-23"),
    Date("2021-12-29"), Date("2021-12-30"), Date("2021-12-31"),
]

const predict_date_start = confirmed_df.date_mid[end]
const predict_date_end   = predict_date_start + Day(60)
@info predict_date_start, predict_date_end
function predict_confirmed(confirmed_df::DataFrame; β)
    @info "--------"
    @info "データフレーム作成"
    df = DataFrame()
    n_rows = Dates.value(predict_date_end - predict_date_start) + 1
    @info "n_rows", n_rows
    @info "--------"
    @info "日付 (date, day_of_week, date_value)"
    df.date = predict_date_start:Day(1):predict_date_end
    df.date_value = date_to_value.(df.date)
    df.day_of_week = dayofweek.(df.date)
    holidays = Vector{Bool}(undef, n_rows + 1)
    for i in 1:n_rows+1
        if i ≤ 1
            d = df.date[1] - Day(1)
            w = dayofweek(d)
        else
            d = df.date[i-1]
            w = df.day_of_week[i-1]
        end
        holidays[i] = (d in holiday_list) && (w != 6) && (w != 7)
    end
    df.holiday = holidays[2:end]
    df.holiday_next = holidays[1:end-1]
    @info size(df)
    @info "--------"
    @info "δ 株/非 δ 株拡大係数設定 (lgrowth_δ_lagXX, lgrowth_nonδ_lagXX)"
    @info "過去 7 日間の平均. 値は時間に依らない定数"
    c_n_rows = nrow(confirmed_df)
    period = min(c_n_rows, 7)
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        ξdname = "lgrowth_δ_"*lagname
        ξds_δ = confirmed_df[!,ξdname]
        ξd_δ_f = mean(ξds_δ[(c_n_rows-period+1):c_n_rows])
        @info ξdname, ξd_δ_f
        df[!,ξdname] = fill(ξd_δ_f, n_rows)
    end
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        ξdname = "lgrowth_nonδ_"*lagname
        ξds_nonδ = confirmed_df[!,ξdname]
        ξd_nonδ_f = mean(ξds_nonδ[(c_n_rows-period+1):c_n_rows])
        @info ξdname, ξd_nonδ_f
        df[!,ξdname] = fill(ξd_nonδ_f, n_rows)
    end
    @info size(df)
    @info "--------"
    @info "拡大係数予測 (lgrowth_lagXX)"
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        ts = df.date_value
        ps = δ_prop_maxL.(ts .+ lag)
        ξds_nonδ = df[!,"lgrowth_nonδ_"*lagname]
        ξds = β .* ps .+ ξds_nonδ
        df[!,"lgrowth_"*lagname] = ξds
    end
    @info size(df)
    @info "--------"
    @info "非 δ 株感染確認数予測 (trend_nonδ_lagXX)"
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        ts = df.date_value
        ξds_nonδ = df[!,"lgrowth_nonδ_"*lagname]
        x_nonδ_1 = confirmed_df[end,"trend_nonδ_"*lagname]
        a_nonδ_1 = x_nonδ_1 / exp(ξds_nonδ[1] * ts[1])
        xs_nonδ = a_nonδ_1 .* exp.(ξds_nonδ .* ts)
        df[!,"trend_nonδ_"*lagname] = xs_nonδ
    end
    @info size(df)
    @info "--------"
    @info "感染確認数予測 (trend_lagXX)"
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        ts = df.date_value
        ps = δ_prop_maxL.(ts .+ lag)
        xs_nonδ = df[!,"trend_nonδ_"*lagname]
        xs = xs_nonδ ./ (1.0 .- ps)
        df[!,"trend_δ_"*lagname] = xs .- xs_nonδ
        df[!,"trend_"*lagname] = xs
    end
    @info "--------"
    @info "週変動付き予測 (trend_w_lagXX, trend_w_cil_lagXX, trend_w_ciu_lagXX)"
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        tw     = zeros(Float64, n_rows)
        tw_cil = zeros(Float64, n_rows)
        tw_ciu = zeros(Float64, n_rows)
        for i in 1:n_rows
            lx = log(df[i,"trend_"*lagname])
            m::Float64 = 0.0
            s::Float64 = 0.0
            if df.holiday_next[i]
                m = ha_ldec_mean
                s = ha_ldec_sd
            else
                w = df.day_of_week[i]
                m = lsdow_df.mean[w]
                s = lsdow_df.sd[w]
            end
            tw[i]     = exp(lx + m)
            tw_cil[i] = exp(lx + m - ci95s * s)
            tw_ciu[i] = exp(lx + m + ci95s * s)
        end
        df[!,"trend_w_"*lagname]    = tw
        df[!,"trend_w_cil_"*lagname] = tw_cil
        df[!,"trend_w_ciu_"*lagname] = tw_ciu
    end
    @info size(df)
    return df
end

predict_df = predict_confirmed(confirmed_df; β=likelihood_max_β)

@info "========"
@info "プロット"

function relx(pl, r)
    l = Plots.xlims(pl)
    return (1.0-r) * l[1] + r * l[2]
end
function rely(pl, r; logscale=false)
    l = Plots.ylims(pl)
    if logscale
        ll = log.(l)
        return exp((1.0-r) * ll[1] + r * ll[2])
    else
        return (1.0-r) * l[1] + r * l[2]
    end
end
#
pyplot(
    titlefont=font("Meiryo",10),
    guidefont=font("Meiryo",10),
    tickfont=font("Meiryo",9)
)
const plot_date_start = confirmed_df.date_mid[end] - Day(60)
const plot_date_end   = plot_date_start + Day(100)
#const plot_date_start = confirmed_df.date_mid[end] - Day(50) - Day(30)
#const plot_date_end   = plot_date_start + Day(100) + Day(30)
const date_ticks_dict = Dict(
    Date("2021-03-15") => "03\u201015",
    Date("2021-03-21") => "緊急事態宣言解除",
    Date("2021-04-01") => "2021\u201004\u201001",
    Date("2021-04-15") => "04\u201015",
    Date("2021-04-25") => "緊急事態宣言発令",
    Date("2021-05-01") => "2021\u201005\u201001",
    Date("2021-05-15") => "05\u201015",
    Date("2021-06-01") => "2021\u201006\u201001",
    Date("2021-06-15") => "06\u201015",
    Date("2021-06-20") => "緊急事態宣言解除",
    Date("2021-07-01") => "2021\u201007\u201001",
    Date("2021-07-04") => "都議会議員選挙",
    Date("2021-07-12") => "緊急事態宣言発令",
    Date("2021-07-15") => "07\u201015",
    Date("2021-07-20") => "1\\ 学期終業式",
    Date("2021-07-23") => "オリンピック開会",
    Date("2021-08-01") => "2021\u201008\u201001",
    Date("2021-08-08") => "オリンピック閉会",
    Date("2021-08-15") => "08\u201015",
    Date("2021-09-01") => "2021\u201009\u201001",
)
const date_ticks_dates = sort(collect(keys(date_ticks_dict)))
const date_ticks_dict_en = Dict(
    Date("2021-03-21") => "End\\ of\\ the\\ SoE",
    Date("2021-04-01") => "2021\u201004\u201001",
    Date("2021-04-15") => "04\u201015",
    Date("2021-04-25") => "Decl.\\ of\\ the\\ SoE",
    Date("2021-05-01") => "2021\u201005\u201001",
    Date("2021-05-15") => "05\u201015",
    Date("2021-06-01") => "2021\u201006\u201001",
    Date("2021-06-15") => "06\u201015",
    Date("2021-06-20") => "End\\ of\\ the\\ SoE",
    Date("2021-07-01") => "2021\u201007\u201001",
    Date("2021-07-04") => "Assembly\\ election",
    Date("2021-07-12") => "Decl.\\ of\\ the\\ SoE",
    Date("2021-07-15") => "07\u201015",
    Date("2021-07-20") => "End\\ of\\ 1st\\ sem.",
    Date("2021-07-23") => "Olympics\\ starts",
    Date("2021-08-01") => "2021\u201008\u201001",
    Date("2021-08-08") => "Olympics\\ ends",
    Date("2021-08-15") => "08\u201015",
    Date("2021-09-01") => "2021\u201009\u201001",
)
const date_ticks_dates_en = sort(collect(keys(date_ticks_dict_en)))

plots = Dict{String,Plots.Plot}()

@info "--------"
@info "確認感染者数・トレンド"
for scale in ["normal", "log"]
    for lang in ["ja", "en"]
        for lag in likelihood_lags
            lagname = @sprintf("lag%02d", lag)
            title = (lang=="ja") ?
                "東京（変異株報告における追加遅延 $(lag) 日）" :
                "Tōkyō (additional $(lag)-day delay in variant reporting)"
            yscale = (scale=="normal") ? :none : :log
            pl = plot(
                size=(640, 700),
                framestyle=:box,
                bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
                title=title,
                legend=:topleft,
                fontfamily="Meiryo",
                yscale=yscale
            )
            # x 軸設定
            sdate_value = date_to_value(plot_date_start)
            edate_value = date_to_value(plot_date_end)
            xlims!(pl, sdate_value, edate_value)
            nticks = length(date_ticks_dates)
            vdates = Float64[]
            labels = String[]
            tdict = (lang=="ja") ? date_ticks_dict : date_ticks_dict_en
            for i in 1:nticks
                date = date_ticks_dates[i]
                if plot_date_start ≤ date ≤ plot_date_end
                    push!(vdates, date_to_value(date))
                    push!(labels, tdict[date])
                end
            end
            xticks!(pl, vdates, labels; xtickfontrotation=90.0)
            # y 軸設定
            if scale == "normal"
                ylims!(pl, 0.0, 5000.0),
                ylabel!(pl, (lang=="ja") ? "感染確認者数 [人/日]" : "Confirmed cases [cases/day]")
                yticks!(pl,
                    [0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0],
                    ["0", "1000", "2000", "3000", "4000", "5000"]
                )
            else
                ylims!(pl, 150.0, 10000.0),
                ylabel!(pl, (lang=="ja") ? "感染確認者数 [人/日]（対数目盛り）" : "Confirmed cases [cases/day] (logscale)")
                yticks!(pl,
                    [150.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0],
                    ["150", "200", "400", "600", "800", "1000", "1500", "2000", "4000", "6000", "8000", "10000"]
                )
            end
            # 感染確認者数
            ts_c   = confirmed_df.date_value
            cs_c   = confirmed_df.count
            xs_c   = confirmed_df.trend
            xs_δ_c = confirmed_df[!,"trend_δ_"*lagname]
            xs_o_c = confirmed_df[!,"trend_nonδ_"*lagname]
            plot!(pl, ts_c, cs_c;
                label= (lang=="ja") ? "感染確認者数" : "confirmed cases",
                lc=:royalblue,
                m=:circle, mc=:royalblue1
            )
            plot!(pl, ts_c, xs_c;
                label= (lang=="ja") ? "トレンド（LOESS 近似）" : "trend (LOESS approx.)",
                lc=:red3, lw=2
            )
            plot!(pl, ts_c, xs_o_c;
                label= (lang=="ja") ? "非デルタ株（推定・予測）" : "non-δ variant (est./pred.)",
                lc=:cyan,
                frange=zeros(Float64, nrow(confirmed_df)),
                fa=0.3, fc=:cyan
            )
            if scale=="normal"
                plot!(pl, ts_c, xs_c;
                    label= (lang=="ja") ? "デルタ株（推定・予測）" : "δ variant (est./pred.)",
                    lc=:red3,
                    frange=xs_o_c,
                    fa=0.4, fc=:orange
                )
            else
                plot!(pl, ts_c, xs_δ_c;
                    label= (lang=="ja") ? "デルタ株（推定・予測）" : "δ variant (est./pred.)",
                    lc=:orange,
                    frange=zeros(Float64, nrow(confirmed_df)),
                    fa=0.4, fc=:orange
                )
            end
            plot!(pl, ts_c, xs_o_c; label=:none, lc=:cyan)
            if scale=="log"
                plot!(pl, ts_c, xs_δ_c; label=:none, lc=:orange)
            end
            plot!(pl, ts_c, xs_c;   label=:none, lc=:red3, lw=2)
            plot!(pl, ts_c, cs_c;   label=:none, lc=:royalblue, m=:circle, mc=:royalblue1)
            # 感染確認者予測
            ts_p    = predict_df.date_value
            xs_p    = predict_df[!,"trend_"*lagname]
            xs_w_p  = predict_df[!,"trend_w_"*lagname]
            xs_wl_p = predict_df[!,"trend_w_cil_"*lagname]
            xs_wu_p = predict_df[!,"trend_w_ciu_"*lagname]
            xs_δ_p  = predict_df[!,"trend_δ_"*lagname]
            xs_o_p  = predict_df[!,"trend_nonδ_"*lagname]
            if scale == "normal"
                plot!(pl, ts_p, xs_p;
                    label=:none,
                    lc=:red3, ls=:dot, lw=2,
                    frange=xs_o_p,
                    fa=0.1, fc=:orange
                )
            else
                plot!(pl, ts_p, xs_δ_p;
                    label=:none,
                    lc=:orange, ls=:dot,
                    frange=zeros(Float64, nrow(predict_df)),
                    fa=0.1, fc=:orange
                )
            end
            plot!(pl, ts_p, xs_o_p;
                label=:none,
                lc=:cyan, ls=:dot,
                frange=zeros(Float64, nrow(predict_df)),
                fa=0.1, fc=:cyan
            )
            #=
            plot!(pl, ts_p[2:end], xs_wu_p[2:end];
                label=:none,
                la=0.0,
                frange=xs_wl_p[2:end],
                fa=0.1, fc=:royalblue
            )
            =#
            plot!(pl, ts_p[2:end], xs_w_p[2:end];
                label= (lang=="ja") ? "感染確認者公表数予測" : "pred. # of confirmed cases",
                la=0.2, lc=:royalblue, ls=:dot,
                m=:circle, ma=0.2, mc=:royalblue1
            )
            plot!(pl, ts_p, xs_p;
                label=:none,
                lc=:red3, ls=:dot, lw=2
            )
            # アノテーション
            islog = (scale == "log")
            annotate!(pl, relx(pl, 0.98), rely(pl, 0.05; logscale=islog),
                text(
                    (lang=="ja") ?
                    """
                    予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（各株に対し係数一定）と想定
                    デルタ株確認者数の予測は最尤線をたどるとした場合
                    予測値は変異株の報告における追加の遅延日数をどれだけとみるかにより変化する
                    """ :
                    """
                    The predictions are based on the expansion coefficient for the recent week,
                    and assume that the equivalent response would continue.
                    The prediction of the δ cases assumes that the maximum likelihood line is followed.
                    The predictions vary depending on additional delay in reporting the variants.
                    """,
                    font("Meiryo",7), RGB(0.3,0.3,0.3), :right
                )
            )
            plname = "confirmed_"*lagname*"_"*scale*"_"*lang
            plots[plname] = pl
            @info plname
        end
    end
end
@info "--------"
@info "確認感染者数・トレンド（ラグ統合）"
for scale in ["normal", "log"]
    for lang in ["ja", "en"]
        yscale = (scale=="normal") ? :none : :log
        pl = plot(
            size=(640,700),
            framestyle=:box,
            bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
            title= (lang=="ja") ? "東京" : "Tōkyō",
            legend=:topleft,
            fontfamily="Meiryo",
            yscale=yscale
        )
        # x 軸設定
        sdate_value = date_to_value(plot_date_start)
        edate_value = date_to_value(plot_date_end)
        xlims!(pl, sdate_value, edate_value)
        nticks = length(date_ticks_dates)
        vdates = Float64[]
        labels = String[]
        tdict = (lang=="ja") ? date_ticks_dict : date_ticks_dict_en
        for i in 1:nticks
            date = date_ticks_dates[i]
            if plot_date_start ≤ date ≤ plot_date_end
                push!(vdates, date_to_value(date))
                push!(labels, tdict[date])
            end
        end
        xticks!(pl, vdates, labels; xtickfontrotation=90.0)
        # y 軸設定
        if scale == "normal"
            ylims!(pl, 0.0, 5000.0)
            ylabel!(pl, (lang=="ja") ? "感染確認者数 [人/日]" : "Confirmed cases [cases/day]")
            yticks!(pl,
                [0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0],
                ["0", "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000", "4500", "5000"]
            )
        else
            ylims!(pl, 150.0, 8000.0)
            ylabel!(pl, (lang=="ja") ? "感染確認者数 [人/日]（対数目盛り）" : "Confirmed cases [cases/day] (logscale)")
            yticks!(pl,
                [150.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 4000.0, 6000.0, 8000.0],
                ["150", "200", "400", "600", "800", "1000", "1500", "2000", "4000", "6000", "8000"]
            )
        end
        # 感染確認者数
        ts_c = confirmed_df.date_value
        cs_c = confirmed_df.count
        xs_c = confirmed_df.trend
        plot!(pl, ts_c, cs_c;
            label= (lang=="ja") ? "感染確認者数" : "confirmed cases",
            lc=:royalblue,
            m=:circle, mc=:royalblue1
        )
        plot!(pl, ts_c, xs_c;
            label= (lang=="ja") ? "トレンド" : "trend",
            lc=:red3, lw=2
        )
        xs_o_c = confirmed_df.trend_nonδ_lag00
        plot!(pl, ts_c, xs_o_c;
            label= (lang=="ja") ? "非デルタ株（推定・予測）" : "non-δ variant (est./pred.)",
            lc=:cyan,
            frange=zeros(Float64, nrow(confirmed_df)),
            fa=0.2, fc=:cyan
        )
        plot!(pl, ts_c, xs_c;
            label= (lang=="ja") ? "デルタ株（推定・予測）" : "δ variant (est./pred.)",
            lc=:red3,
            frange=xs_o_c,
            fa=0.2, fc=:orange
        )
        plot!(pl, ts_c, xs_c;
            label= (lang=="ja") ? "デルタ株（推定・予測）" : "δ variant (est./pred.)",
            lc=:red3,
        )
        for lag in likelihood_lags
            lagname = @sprintf("lag%02d", lag)
            xs_l_c = confirmed_df[!,"trend_nonδ_"*lagname]
            if lag > 0
                plot!(pl, ts_c, xs_o_c;
                    label=:none,
                    la=0.0,
                    frange=xs_l_c,
                    fa=0.2, fc=:orange
                )
                plot!(pl, ts_c, xs_l_c;
                    label=:none,
                    lc=:cyan,
                    frange=zeros(Float64, nrow(confirmed_df)),
                    fa=0.2, fc=:cyan
                )
            end        
        end
        plot!(pl, ts_c, xs_c; label=:none, lc=:red3, lw=2)
        plot!(pl, ts_c, cs_c; label=:none, lc=:royalblue, m=:circle, mc=:royalblue1)
        # 感染確認者予測
        ts_p   = predict_df.date_value
        xs_p   = predict_df.trend_lag00
        xs_o_p = predict_df.trend_nonδ_lag00
        plot!(pl, ts_p, xs_p;
            label=:none,
            lc=:red3, ls=:dot, lw=2,
            frange=xs_o_p,
            fa=0.05, fc=:orange
        )
        plot!(pl, ts_p, xs_o_p;
            label=:none,
            lc=:cyan, ls=:dot,
            frange=zeros(Float64, nrow(predict_df)),
            fa=0.05, fc=:cyan
        )
        for lag in likelihood_lags
            lagname = @sprintf("lag%02d", lag)
            xs_l_p = predict_df[!,"trend_nonδ_"*lagname]
            if lag > 0
                plot!(pl, ts_p, xs_o_p;
                    label=:none,
                    la=0.0,
                    frange=xs_l_p,
                    fa=0.05, fc=:orange
                )
                plot!(pl, ts_p, xs_l_p;
                    label=:none,
                    la=0.3, lc=:cyan, ls=:dot,
                    frange=zeros(Float64, nrow(predict_df)),
                    fa=0.05, fc=:cyan
                )
            end
        end
        # アノテーション
        islog = (scale == "log")
        annotate!(pl, relx(pl, 0.98), rely(pl, 0.05; logscale=islog),
            text(
                (lang=="ja") ?
                """
                予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
                デルタ株確認者数の予測は最尤線をたどるとした場合
                予測値は変異株の報告における追加の遅延日数をどれだけとみるかにより変化する
                """ :
                """
                The predictions are based on the expansion coefficient for the recent week,
                and assume that the equivalent response would continue.
                The prediction of the δ cases assumes that the maximum likelihood line is followed.
                The predictions vary depending on additional delay in reporting the variants.
                """,
                font("Meiryo",7), RGB(0.3,0.3,0.3), :right
            )
        )
        plname = "confirmed_"*scale*"_"*lang
        plots[plname] = pl
        @info plname
    end
end
@info "--------"
@info "拡大速度係数"
for scale in ["normal", "week", "doubling", "reprod"]
    for lang in ["ja", "en"]
        for lag in likelihood_lags
            lagname = @sprintf("lag%02d", lag)
            if lag == 0
                title = (lang=="ja") ? "東京" : "Tōkyō"
            else
                title = (lang=="ja") ?
                    "東京（変異株報告における追加遅延 $(lag) 日）" :
                    "Tōkyō (additional $(lag)-day delay in variant reporting)"
            end
            pl = plot(
                size=(640, 600),
                framestyle=:box,
                bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
                title=title,
                legend=:topleft,
                fontfamily="Meiryo"
            )
            # x 軸設定
            sdate_value = date_to_value(plot_date_start)
            edate_value = date_to_value(plot_date_end)
            xlims!(pl, sdate_value, edate_value)
            nticks = length(date_ticks_dates)
            vdates = Float64[]
            labels = String[]
            tdict = (lang=="ja") ? date_ticks_dict : date_ticks_dict_en
            for i in 1:nticks
                date = date_ticks_dates[i]
                if plot_date_start ≤ date ≤ plot_date_end
                    push!(vdates, date_to_value(date))
                    push!(labels, tdict[date])
                end
            end
            xticks!(pl, vdates, labels; xtickfontrotation=90.0)
            # y 軸設定
            ymin = -0.1
            ymax =  0.2
            ylims!(pl, ymin, ymax)
            if     scale=="normal"
                ylabel!(pl, (lang=="ja") ? "拡大係数（時定数の逆数）[/日]" : "Expansion coefficient (1 / time constant) [/day]")
                ecs = ymin:0.05:ymax
                yticks!(pl, ecs, [@sprintf("%.2f", ec) for ec in ecs])
            elseif scale=="week"
                ylabel!(pl, (lang=="ja") ? "週あたり拡大率（対数目盛り）" : "Expansion rate per week (logscale)")
                wrmin = 0.2 * ceil(exp(7.0 * ymin) / 0.2)
                wrmax = 0.5 * floor(exp(7.0 * ymax) / 0.5)
                wrs = [wrmin:0.2:1.8; 2.0:0.5:wrmax]
                yticks!(pl,
                    [log(wr)/7.0 for wr in wrs],
                    [@sprintf("%.1f", wr) for wr in wrs]
                )
            elseif scale=="doubling"
                ylabel!(pl, (lang=="ja") ? "倍加時間 [日]" : "Doubling time [day]")
                dt0s = [2.0, 3.0, 4.0, 5.0, 7.0, 14.0, 28.0]
                #       0.50 0.33 0.25 0.20 0.14 0.07  0.04
                dtvs = Float64[]
                dtls = String[]
                for dt0 in dt0s
                    ec = - log(2.0) / dt0
                    if ymin ≤ ec ≤ ymax
                        push!(dtvs, ec)
                        push!(dtls, @sprintf("%.0f", -dt0))
                    end
                end
                push!(dtvs, 0.0)
                push!(dtls, "∞")
                for dt0 in reverse(dt0s)
                    ec = log(2.0) / dt0
                    if ymin ≤ ec ≤ ymax
                        push!(dtvs, ec)
                        push!(dtls, @sprintf("%.0f", dt0))
                    end
                end
                yticks!(pl, dtvs, dtls)
            end
            # グリッド
            plot!(pl, [sdate_value, edate_value], [0.0, 0.0]; label=:none, lc=RGB(0.3,0.3,0.3))
            # 拡大係数
            ts_c    = confirmed_df.date_value[1:end-4]
            ξds_c   = confirmed_df.lgrowth[1:end-4]
            ξds_δ_c = (confirmed_df[!,"lgrowth_δ_"*lagname])[1:end-4]
            ξds_o_c = (confirmed_df[!,"lgrowth_nonδ_"*lagname])[1:end-4]
            if     scale=="normal"
                eclbl = (lang=="ja") ? "拡大係数" : "Expansion coefficient"
            elseif scale=="week"
                eclbl = (lang=="ja") ? "週あたり拡大率" : "Expansion rate per week"
            elseif scale=="doubling"
                eclbl = (lang=="ja") ? "倍加時間" : "Doubling time"
            else
                eclbl = ""
            end
            plot!(pl, ts_c, ξds_c;
                label=eclbl,
                lc=:red3, lw=2
            )
            plot!(pl, ts_c, ξds_δ_c;
                label= (lang=="ja") ? "デルタ株推定" : "δ variant estimation",
                lc=:orange
            )
            plot!(pl, ts_c, ξds_o_c;
                label= (lang=="ja") ? "非デルタ株推定" : "non-δ variant estimation",
                lc=:cyan
            )
            plot!(pl, ts_c, ξds_c; label=:none, lc=:red3, lw=2)
            # 拡大係数予測
            ts_p    = predict_df.date_value
            ξds_p   = predict_df[!,"lgrowth_"*lagname]
            ξds_δ_p = predict_df[!,"lgrowth_δ_"*lagname]
            ξds_o_p = predict_df[!,"lgrowth_nonδ_"*lagname]
            plot!(pl, ts_p, ξds_δ_p;
                label= (lang=="ja") ? "デルタ株想定値" : " δ variant assumption",
                lc=:orange, ls=:dot, lw=2
            )
            plot!(pl, ts_p, ξds_o_p;
                label= (lang=="ja") ? "非デルタ株想定値" : "non-δ variant assumption",
                lc=:cyan, ls=:dot, lw=2
            )
            plot!(pl, ts_p, ξds_p;
                label= (lang=="ja") ? "デルタ株を加えた予測" : "prediction with δ variant",
                lc=:red3, ls=:dot, lw=2
            )
            # アノテーション
            if     scale=="normal"
                anntxt = (lang=="ja") ?
                    """
                    ここで拡大係数とは確認者数に指数関数をあてはめたときの時定数の逆数. 正で拡大、負で減少
                    拡大係数が正のとき、その逆数の約 0.7 倍が拡大の倍加時間にあたる（0.05 なら倍加時間約 14 日）
                    予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
                    デルタ株拡大の予測は最尤線をたどるとした場合
                    """ :
                    """
                    “Expansion coefficient” denotes here the inverse of local time constants of confirmed cases.
                    If the coefficient is positive, the doubling time is ~0.7 times of the inverse of it.
                    The prediction of the δ cases assumes that the maximum likelihood line is followed.
                    The predictions vary depending on additional delay in reporting the variants.
                    """ 
            else
                anntxt = (lang=="ja") ?
                    """
                    予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
                    デルタ株拡大の予測は最尤線をたどるとした場合
                    """ :
                    """
                    The prediction of the δ cases assumes that the maximum likelihood line is followed.
                    The predictions vary depending on additional delay in reporting the variants.
                    """ 
            end
            annotate!(pl, relx(pl, 0.02), rely(pl, 0.05),
                text(anntxt, font("Meiryo",7), RGB(0.3,0.3,0.3), :left)
            )
            plname = "growth_"*lagname*"_"*scale*"_"*lang
            plots[plname] = pl
            @info plname
        end
    end
end

@info "--------"
@info "ゴールデン・ウィーク"
const date_ticks_dict_gw = Dict(
    Date("2021-04-03") => "2021\u201004\u201003（土）",
    Date("2021-04-04") => "2021\u201004\u201004（日）",
    Date("2021-04-10") => "2021\u201004\u201010（土）",
    Date("2021-04-11") => "2021\u201004\u201011（日）",
    Date("2021-04-17") => "2021\u201004\u201017（土）",
    Date("2021-04-18") => "2021\u201004\u201018（日）",
    Date("2021-04-24") => "2021\u201004\u201024（土）",
    Date("2021-04-25") => "2021\u201004\u201025（日）",
    Date("2021-04-29") => "昭和の日（祝）",
    Date("2021-05-01") => "2021\u201005\u201001（土）",
    Date("2021-05-02") => "2021\u201005\u201002（日）",
    Date("2021-05-03") => "憲法記念日（祝）",
    Date("2021-05-04") => "みどりの日（祝）",
    Date("2021-05-05") => "こどもの日（祝）",
    Date("2021-05-08") => "2021\u201005\u201008（土）",
    Date("2021-05-09") => "2021\u201005\u201009（日）",
    Date("2021-05-15") => "2021\u201005\u201015（土）",
    Date("2021-05-16") => "2021\u201005\u201016（日）",
    Date("2021-05-22") => "2021\u201005\u201022（土）",
    Date("2021-05-23") => "2021\u201005\u201023（日）",
    Date("2021-05-29") => "2021\u201005\u201029（土）",
    Date("2021-05-30") => "2021\u201005\u201030（日）",
)
const date_ticks_dates_gw = sort(collect(keys(date_ticks_dict_gw)))
for once in 1:1
    pl = plot(
        size=(720,600),
        framestyle=:box,
        bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
        title="東京（ゴールデン・ウィーク連休時の休日効果）",
        legend=:topleft,
        fontfamily="Meiryo"
    )
    # x 軸設定
    sdate_value = date_to_value(Date("2021-04-14"))
    edate_value = date_to_value(Date("2021-05-27"))
    xlims!(pl, sdate_value, edate_value)
    nticks = length(date_ticks_dates_gw)
    vdates = Float64[]
    labels = String[]
    for i in 1:nticks
        date = date_ticks_dates_gw[i]
        if sdate_value ≤ date_to_value(date) ≤ edate_value
            push!(vdates, date_to_value(date) + 0.5)
            push!(labels, date_ticks_dict_gw[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, 0.0, 1400.0)
    ylabel!(pl, "感染確認者数 [人/日]")
    yticks!(pl,
        [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0],
        ["0", "200", "400", "600", "800", "1000", "1200", "1400"]
    )
    # 感染確認者数
    ts_c = confirmed_df.date_value
    cs_c = confirmed_df.count
    xs_c = confirmed_df.trend
    plot!(pl, ts_c, cs_c; label="感染確認者数", m=:circle, lc=:royalblue, mc=:royalblue1)
    plot!(pl, ts_c, xs_c; label="トレンド", lw=2, lc=:red3)
    xs_nonδ_c = confirmed_df.trend_nonδ_lag00
    plot!(pl, ts_c, cs_c; label=:none, m=:circle, lc=:royalblue, mc=:royalblue1)
    plname = "gw_holidays"
    plots[plname] = pl
    @info plname
end
#=
@info "--------"
@info "指数関数の足し合わせ"
for once in 1:1
    pl = plot(
        size=(640, 700),
        framestyle=:box,
        bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
        legend=:topleft,
        fontfamily="Meiryo"
    )
    # x 軸設定
    sdate_value = date_to_value(plot_date_start + Day(20))
    edate_value = date_to_value(plot_date_end - Day(20))
    xlims!(pl, sdate_value, edate_value)
    nticks = length(date_ticks_dates)
    vdates = Float64[]
    labels = String[]
    tdict = date_ticks_dict
    for i in 1:nticks
        date = date_ticks_dates[i]
        if sdate_value ≤ date_to_value(date) ≤ edate_value
            push!(vdates, date_to_value(date))
            push!(labels, tdict[date])
        end
        push!(vdates, 198.0)
        push!(labels, "t = 0")
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, 0.0, 4000.0),
    ylabel!(pl, "感染確認者数 [人/日]")
    yticks!(pl,
        [0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0],
        ["0", "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000"]
    )
    # 感染確認者数
    ts_c   = confirmed_df.date_value
    cs_c   = confirmed_df.count
    ts   = sdate_value:1.0:edate_value
    xs_δ = [600.0*(2)^((t-198.0)/7.0) for t in ts]
    xs_o = [600.0*(1.2)^((t-198.0)/7.0) for t in ts]
    xs   = xs_δ .+ xs_o
    plot!(pl, ts_c, cs_c;
        label="感染確認者数",
        lc=:royalblue,
        m=:circle, mc=:royalblue1
    )
    plot!(pl, ts, xs_o;
        label="α(t) = 600・1.2^(t/7)",
        lc=:cyan,
        frange=zeros(Float64, nrow(confirmed_df)),
        fa=0.3, fc=:cyan
    )
    plot!(pl, ts, xs_δ;
        label="δ(t) = 600・2^(t/7)",
        lc=:orange,
        frange=zeros(Float64, nrow(confirmed_df)),
        fa=0.4, fc=:orange
    )
    plot!(pl, ts, xs;
        label="α(t) + δ(t)",
        lc=:red3, lw=2
    )
    plot!(pl, ts, xs_o; label=:none, lc=:cyan)
    plot!(pl, ts, xs_δ; label=:none, lc=:orange)
    plot!(pl, ts, xs;   label=:none, lc=:red3, lw=2)
    plot!(pl, ts_c, cs_c;   label=:none, lc=:royalblue, m=:circle, mc=:royalblue1)
    plname = "math"
    plots[plname] = pl
    @info plname
end
=#
@info "--------"
@info "だめ予測"
for scale in ["log"]
    for lang in ["ja"]
        for lag in [7]
            lagname = @sprintf("lag%02d", lag)
            title = (lang=="ja") ?
                "東京（変異株報告における追加遅延 $(lag) 日）" :
                "Tōkyō (additional $(lag)-day delay in variant reporting)"
            yscale = (scale=="normal") ? :none : :log
            pl = plot(
                size=(640, 700),
                framestyle=:box,
                bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
                title=title,
                legend=:topleft,
                fontfamily="Meiryo",
                yscale=yscale
            )
            # x 軸設定
            sdate_value = date_to_value(plot_date_start)
            edate_value = date_to_value(plot_date_end)
            xlims!(pl, sdate_value, edate_value)
            nticks = length(date_ticks_dates)
            vdates = Float64[]
            labels = String[]
            tdict = (lang=="ja") ? date_ticks_dict : date_ticks_dict_en
            for i in 1:nticks
                date = date_ticks_dates[i]
                if plot_date_start ≤ date ≤ plot_date_end
                    push!(vdates, date_to_value(date))
                    push!(labels, tdict[date])
                end
            end
            xticks!(pl, vdates, labels; xtickfontrotation=90.0)
            # y 軸設定
            if scale == "normal"
                ylims!(pl, 0.0, 4000.0),
                ylabel!(pl, (lang=="ja") ? "感染確認者数 [人/日]" : "Confirmed cases [cases/day]")
                yticks!(pl,
                    [0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0],
                    ["0", "500", "1000", "1500", "2000", "2500", "3000", "3500", "4000"]
                )
            else
                ylims!(pl, 150.0, 12000.0),
                ylabel!(pl, (lang=="ja") ? "感染確認者数 [人/日]（対数目盛り）" : "Confirmed cases [cases/day] (logscale)")
                yticks!(pl,
                    [150.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 4000.0, 6000.0, 8000.0, 10000.0],
                    ["150", "200", "400", "600", "800", "1000", "1500", "2000", "4000", "6000", "8000", "10000"]
                )
            end
            # 感染確認者数
            ts_c   = confirmed_df.date_value
            cs_c   = confirmed_df.count
            xs_c   = confirmed_df.trend
            xs_δ_c = confirmed_df[!,"trend_δ_"*lagname]
            xs_o_c = confirmed_df[!,"trend_nonδ_"*lagname]
            plot!(pl, ts_c, cs_c;
                label= (lang=="ja") ? "感染確認者数" : "confirmed cases",
                la=0.5, lc=:royalblue,
                m=:circle, ma=0.7, mc=:royalblue1
            )
            plot!(pl, ts_c, xs_c;
                label= (lang=="ja") ? "トレンド（LOESS 近似）" : "trend (LOESS approx.)",
                lc=:red3, lw=2
            )
            plot!(pl, ts_c, xs_c; label=:none, lc=:red3, lw=2)
            # 感染確認者予測
            ts_p    = predict_df.date_value
            xs_p    = predict_df[!,"trend_"*lagname]
            xs_w_p  = predict_df[!,"trend_w_"*lagname]
            xs_wl_p = predict_df[!,"trend_w_cil_"*lagname]
            xs_wu_p = predict_df[!,"trend_w_ciu_"*lagname]
            xs_δ_p  = predict_df[!,"trend_δ_"*lagname]
            xs_o_p  = predict_df[!,"trend_nonδ_"*lagname]
            plot!(pl, ts_p, xs_p;
                label=:none,
                lc=:red3, ls=:dot, lw=2
            )
            # 都モニタリング会議によるだめ予測
            todame_preds = [
                [(Date("2021-06-30"),  503), (Date("2021-06-30")+Day(4*7), 1043)],
                [(Date("2021-07-07"),  625), (Date("2021-07-07")+Day(4*7), 1478)],
                [(Date("2021-07-14"),  817), (Date("2021-07-14")+Day(4*7), 2406)],
                [(Date("2021-07-20"), 1170), (Date("2021-07-20")+Day(2*7), 2598)],                
                [(Date("2021-07-28"), 1936), (Date("2021-07-28")+Day(2*7), 4582)],                
            ]
            tdm_label::Union{String,Symbol} = "都モニタリング会議試算"
            for dps in todame_preds
                plot!(pl,
                    [date_to_value(dp[1]) for dp in dps],
                    [dp[2] for dp in dps],
                    label=tdm_label,
                    lc=:green, lw=1,
                    m=:diamond, mc=:forestgreen
                )
                tdm_label=:none
            end
            # 西浦氏によるだめ予測
            nidame_preds = [
                [(Date("2021-08-11"), 5027), (Date("2021-08-26"), 10643)],
            ]
            for dps in nidame_preds
                plot!(pl,
                    [date_to_value(dp[1]) for dp in dps],
                    [dp[2] for dp in dps],
                    label="2021-07-28 アドバイザリーボード 西浦試算",
                    lc=:cyan, ls=:dot, lw=1,
                    m=:diamond, mc=:cyan
                )
            end
            # アノテーション
            islog = (scale == "log")
            annotate!(pl, relx(pl, 0.98), rely(pl, 0.05; logscale=islog),
                text(
                    (lang=="ja") ?
                    """
                    予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（各株に対し係数一定）と想定
                    デルタ株確認者数の予測は最尤線をたどるとした場合
                    予測値は変異株の報告における追加の遅延日数をどれだけとみるかにより変化する
                    """ :
                    """
                    The predictions are based on the expansion coefficient for the recent week,
                    and assume that the equivalent response would continue.
                    The prediction of the δ cases assumes that the maximum likelihood line is followed.
                    The predictions vary depending on additional delay in reporting the variants.
                    """,
                    font("Meiryo",7), RGB(0.3,0.3,0.3), :right
                )
            )
            plname = "monpred"
            plots[plname] = pl
            @info plname
        end
    end
end

end #module
