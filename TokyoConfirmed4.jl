module TokyoConfirmed4
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
    df.holiday_next = ([false; df.holiday])[1:end-1]
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
const qc_period = 21
@info "qc_period", qc_period

linreg(x, β) = β[1] .+ β[2] .* x
linreg_inv(y, β) = (y .- β[1]) ./ β[2]
qc_reg(x, α) = α[1] .+ x .* (α[2] .+ x .* α[3])

function ssanalysis!(df::DataFrame)
    n_rows = nrow(df)
    @info "--------"
    @info "予備特異スペクトル解析"
    @info "トレンド (xs0), 周期成分 (ss0), 残差 (rs0) への分解. ps0 は ss0 + rs0"
    cs = df.count
    #=
    (xs0, ss0_array) = analyze(cs, ssa_window; robust=false)
    @info "****", size(ss0_array)
    ps0 = cs .- xs0     # pecular = count - trend
    ss0 = zeros(Float64, n_rows)
    for ss0_col in eachcol(ss0_array)
        ss0 .+= ss0_col
    end
    rs0 = ps0 .- ss0    # residual = pecular - seasonal
    =#
    cs_stl = rcopy(R"""
    cs.ts <- ts($cs, frequency=$ssa_window, start=1)
    cs.stl <- stl(cs.ts, s.window=$ssa_window)
    """)
    cs_rts = cs_stl[:time_series]
    @info size(cs_rts)
    ps0 = cs_rts[:,1] .+ cs_rts[:,3]    # seasonal + residual
    @info "--------"
    @info "超過量補正 (ladjust1)"
    @info "周期的変動の大きさが count に依存することの補正"
    ps0_abs = 2.0 * abs.(ps0)
    β0 = [50.0, 0.1]
    p1fit = curve_fit(linreg, cs, ps0_abs, β0)
    @info p1fit.param
    @assert p1fit.param[1] ≥ 0.0
    lcs1 = log.(linreg(cs, p1fit.param))
    df.ladjust1 = lcs1
    @info size(df)
    #=
    @info "--------"
    @info "2次係数補正 (ladjust2)"
    @info "直近の2次係数が0となるよう補正" 
    α0 = [0.0, 0.0, 0.0]
    qc_x = df.date_value[end-qc_period:end]
    qc_y = lcs1[end-qc_period:end]
    p2fit = curve_fit(qc_reg, qc_x, qc_y, α0)
    qc_α2 = p2fit.param[3]
    lcs2 = lcs1 .- qc_α2 .* (df.date_value).^2
    @info "qc_α2", qc_α2
    df.ladjust2 = lcs2
    @info size(df)
    =#
    @info size(df)
    @info "--------"
    @info "特異スペクトル解析 (ltrend1, lseasonal)"
    #=
    (lxs, lss_array) = analyze(lcs2, ssa_window; robust=false)
    @info "****", size(lss_array)
    lps = lcs2 .- lxs
    lss = zeros(Float64, n_rows)
    for lss_col in eachcol(lss_array)
        lss .+= lss_col
    end
    lrs = lps .- lss
    df.ltrend    = lxs
    df.lpeculiar = lps
    df.lseasonal = lss
    df.lresidual = lrs
    =#
    lcs_stl = rcopy(R"""
    lcs.ts <- ts($(df.ladjust1), frequency=$ssa_window, start=1)
    lcs.stl <- stl(lcs.ts, s.window=$ssa_window, l.degree=1, t.degree=1, s.degree=1, robust=T)
    plot(lcs.stl)
    lcs.stl
    """)
    lcs_rts = lcs_stl[:time_series]
    @info size(lcs_rts)
    df.lseasonal = lcs_rts[:,1]
    df.ltrend1   = lcs_rts[:,2]
    df.lresidual = lcs_rts[:,3]
    @info size(df)
    @info "--------"
    @info "残差誤差 (ltrend1_cil, ltrend1_ciu)"
    lrs = filter((row) -> (!row.holiday_next), df).lresidual
    lrm = mean(lrs)
    lre = ci95s * stdm(lrs, lrm)
    @info lrm, lre
    df.ltrend1_cil = df.ltrend1 .- (lre + lrm)
    df.ltrend1_ciu = df.ltrend1 .+ (lre - lrm)
    @info size(df)
    @info "--------"
    @info "補正トレンド (trend, trend_cil, trend_ciu)"
    #=
    df.ltrend1     = df.ltrend2     .+ qc_α2 .* (df.date_value).^2
    df.ltrend1_cil = df.ltrend2_cil .+ qc_α2 .* (df.date_value).^2
    df.ltrend1_ciu = df.ltrend2_ciu .+ qc_α2 .* (df.date_value).^2
    =#
    df.trend     = linreg_inv(exp.(df.ltrend1),     p1fit.param)
    df.trend_cil = linreg_inv(exp.(df.ltrend1_cil), p1fit.param)
    df.trend_ciu = linreg_inv(exp.(df.ltrend1_ciu), p1fit.param)
    for i in 1:n_rows
        @assert df.trend[i] ≥ 0.0
    end
    @info size(df)
    @info "--------"
    @info "拡大係数 (lgrowth, lgrowth_cil, lgrowth_ciu)"
    xs     = df.trend
    xs_cil = df.trend_cil
    xs_ciu = df.trend_ciu
    ξs     = log.(xs)
    ξs_cil = log.(xs_cil)
    ξs_ciu = log.(xs_ciu)
    ξds     = zeros(Float64, n_rows)
    ξds_cil = zeros(Float64, n_rows)
    ξds_ciu = zeros(Float64, n_rows)
    for i in 1:n_rows
        if     i == 1
            ξd = ξs[i+1] - ξs[i]
            ξd_cil = ξd - sqrt((ξs[i+1] - ξs_cil[i+1])^2 + (ξs_ciu[i] - ξs[i])^2)
            ξd_ciu = ξd + sqrt((ξs_ciu[i+1] - ξs[i+1])^2 + (ξs[i] - ξs_cil[i])^2)
        elseif i == n_rows
            ξd = ξs[i] - ξs[i-1]
            ξd_cil = ξd - sqrt((ξs[i] - ξs_cil[i])^2 + (ξs_ciu[i-1] - ξs[i-1])^2)
            ξd_ciu = ξd + sqrt((ξs_ciu[i] - ξs[i])^2 + (ξs[i-1] - ξs_cil[i-1])^2)
        else
            ξd = (ξs[i+1] - ξs[i-1]) / 2.0
            ξd_cil = ξd - sqrt((ξs[i+1] - ξs_cil[i+1])^2 + (ξs_ciu[i-1] - ξs[i-1])^2) / 2.0
            ξd_ciu = ξd + sqrt((ξs_ciu[i+1] - ξs[i+1])^2 + (ξs[i-1] - ξs_cil[i-1])^2) / 2.0
        end
        ξds[i]     = ξd
        ξds_cil[i] = ξd_cil
        ξds_ciu[i] = ξd_ciu
    end
    df.lgrowth     = ξds
    df.lgrowth_cil = ξds_cil
    df.lgrowth_ciu = ξds_ciu
    @info size(df)
end

ssanalysis!(confirmed_df)

@info "========"
@info "曜日変動"
function day_of_week_var(confirmed_df::DataFrame)
    wdf = DataFrame()
    lrs = skipmissing([((df.holiday_next[i]) ? missing : df.lresidual[i]) for i in 1:n_rows])
end

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
    @info "日付 (date, date_value)"
    df.date = predict_date_start:Day(1):predict_date_end
    df.date_value = date_to_value.(df.date)
    @info size(df)
    @info "--------"
    @info "非 δ 株拡大係数設定 (lgrowth_nonδ_lagXX)"
    @info "過去 7 日間の平均. 値は時間に依らない定数"
    c_n_rows = nrow(confirmed_df)
    period = min(c_n_rows, 7)
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
        df[!,"trend_"*lagname] = xs
    end
    @info size(df)
    return df
end

predict_df = predict_confirmed(confirmed_df; β=likelihood_max_β)

@info "========"
@info "プロット"

pyplot(
    titlefont=font("Meiryo",10),
    guidefont=font("Meiryo",10),
    tickfont=font("Meiryo",9)
)
const plot_date_start = confirmed_df.date_mid[end] - Day(50)
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
    Date("2021-07-23") => "オリンピック開会",
    Date("2021-08-01") => "2021\u201008\u201001",
    Date("2021-08-08") => "オリンピック閉会",
    Date("2021-08-15") => "08\u201015",
    Date("2021-08-22") => "宣言解除予定",
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
    Date("2021-07-23") => "Olympics\\ starts",
    Date("2021-08-01") => "2021\u201008\u201001",
    Date("2021-08-08") => "Olympics\\ ends",
    Date("2021-08-15") => "08\u201015",
    Date("2021-08-22") => "SoE\\ ends\\ (sched)",
    Date("2021-09-01") => "2021\u201009\u201001",
)
const date_ticks_dates_en = sort(collect(keys(date_ticks_dict_en)))

plots = Dict{String,Plots.Plot}()
@info "--------"
@info "確認感染者数・トレンド"
for lag in likelihood_lags
    lagname = @sprintf("lag%02d", lag)
    title = "東京"
    if lag != 0
        title *= "（時間遅れ $(lag) 日）"
    end
    pl = plot(
        size=(640,600),
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
    for i in 1:nticks
        date = date_ticks_dates[i]
        if plot_date_start ≤ date ≤ plot_date_end
            push!(vdates, date_to_value(date))
            push!(labels, date_ticks_dict[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, 0.0, 1500.0)
    ylabel!(pl, "感染確認者数 [人/日]")
    # 感染確認者数
    ts_c = confirmed_df.date_value
    cs_c = confirmed_df.count
    xs_c = confirmed_df.trend
    xs_nonδ_c = confirmed_df[!,"trend_nonδ_"*lagname]
    plot!(pl, ts_c, cs_c; label="感染確認者数", m=:circle, lc=:royalblue, mc=:royalblue1)
    plot!(pl, ts_c, xs_c; label="トレンド", lw=2, lc=:red3)
    plot!(pl, ts_c, xs_nonδ_c; label="非デルタ株（推定・予測）",
        frange=zeros(Float64, nrow(confirmed_df)), fa=0.3, c=:cyan)
    plot!(pl, ts_c, xs_c; label="デルタ株（推定・予測）",
        frange=xs_nonδ_c, fa=0.3, c=:red3)
    plot!(pl, ts_c, xs_nonδ_c; label=:none, lc=:cyan)
    plot!(pl, ts_c, xs_c; label=:none, lw=2, lc=:red3)
    plot!(pl, ts_c, cs_c; label=:none, m=:circle, lc=:royalblue, mc=:royalblue1)
    # 感染確認者予測
    ts_p = predict_df.date_value
    xs_p = predict_df[!,"trend_"*lagname]
    xs_nonδ_p = predict_df[!,"trend_nonδ_"*lagname]
    plot!(pl, ts_p, xs_nonδ_p; label=:none,
        frange=zeros(Float64, nrow(predict_df)), fa=0.1, c=:cyan, ls=:dot)
    plot!(pl, ts_p, xs_p; label=:none,
        frange=xs_nonδ_p, fa=0.1, c=:red3, lw=2, ls=:dot)
    plot!(pl, ts_p, xs_nonδ_p; label=:none, lc=:cyan, ls=:dot)
    plot!(pl, ts_p, xs_p; label=:none, lw=2, lc=:red3, ls=:dot)
    # アノテーション
    annotate!(pl, date_to_value(plot_date_start) + 2, 100, text("""
        予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
        デルタ株確認者数の予測は最尤線をたどるとした場合""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left))
    plname = "confirmed_"*lagname
    plots[plname] = pl
    @info plname
end
@info "--------"
@info "in English"
for lag in likelihood_lags
    lagname = @sprintf("lag%02d", lag)
    title = "Tokyo"
    if lag != 0
        title *= "(time lag: $(lag) days)"
    end
    pl = plot(
        size=(640,600),
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
    nticks = length(date_ticks_dates_en)
    vdates = Float64[]
    labels = String[]
    for i in 1:nticks
        date = date_ticks_dates[i]
        if plot_date_start ≤ date ≤ plot_date_end
            push!(vdates, date_to_value(date))
            push!(labels, date_ticks_dict_en[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, 0.0, 1500.0)
    ylabel!(pl, "Confirmed cases [cases/day]")
    # 感染確認者数
    ts_c = confirmed_df.date_value
    cs_c = confirmed_df.count
    xs_c = confirmed_df.trend
    xs_nonδ_c = confirmed_df[!,"trend_nonδ_"*lagname]
    plot!(pl, ts_c, cs_c; label="Confirmed cases", m=:circle, lc=:royalblue, mc=:royalblue1)
    plot!(pl, ts_c, xs_c; label="Trend", lw=2, lc=:red3)
    plot!(pl, ts_c, xs_nonδ_c; label="non-δ cases (est./pred.)",
        frange=zeros(Float64, nrow(confirmed_df)), fa=0.3, c=:cyan)
    plot!(pl, ts_c, xs_c; label="δ cases (est./pred.)",
        frange=xs_nonδ_c, fa=0.3, c=:red3)
    plot!(pl, ts_c, xs_nonδ_c; label=:none, lc=:cyan)
    plot!(pl, ts_c, xs_c; label=:none, lw=2, lc=:red3)
    plot!(pl, ts_c, cs_c; label=:none, m=:circle, lc=:royalblue, mc=:royalblue1)
    # 感染確認者予測
    ts_p = predict_df.date_value
    xs_p = predict_df[!,"trend_"*lagname]
    xs_nonδ_p = predict_df[!,"trend_nonδ_"*lagname]
    plot!(pl, ts_p, xs_nonδ_p; label=:none,
        frange=zeros(Float64, nrow(predict_df)), fa=0.1, c=:cyan, ls=:dot)
    plot!(pl, ts_p, xs_p; label=:none,
        frange=xs_nonδ_p, fa=0.1, c=:red3, lw=2, ls=:dot)
    plot!(pl, ts_p, xs_nonδ_p; label=:none, lc=:cyan, ls=:dot)
    plot!(pl, ts_p, xs_p; label=:none, lw=2, lc=:red3, ls=:dot)
    # アノテーション
    annotate!(pl, date_to_value(plot_date_start) + 2, 100, text("""
        The predictions are based on the expansion coefficient for the recent week
        and assume that the equivalent response would continue.
        The prediction of the δ cases assumes that the maximum likelihood line is followed.""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left))
    plname = "confirmed_"*lagname*"_en"
    plots[plname] = pl
    @info plname
end
@info "--------"
@info "確認感染者数・トレンド（ラグ統合）"
for once in 1:1
    pl = plot(
        size=(640,700),
        framestyle=:box,
        bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
        title="東京",
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
    for i in 1:nticks
        date = date_ticks_dates[i]
        if plot_date_start ≤ date ≤ plot_date_end
            push!(vdates, date_to_value(date))
            push!(labels, date_ticks_dict[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, 0.0, 3500.0)
    ylabel!(pl, "感染確認者数 [人/日]")
    yticks!(pl,
        [0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0],
        ["0", "500", "1000", "1500", "2000", "2500", "3000", "3500"]
    )
    # 感染確認者数
    ts_c = confirmed_df.date_value
    cs_c = confirmed_df.count
    xs_c = confirmed_df.trend
    plot!(pl, ts_c, cs_c; label="感染確認者数", m=:circle, lc=:royalblue, mc=:royalblue1)
    plot!(pl, ts_c, xs_c; label="トレンド", lw=2, lc=:red3)
    xs_nonδ_c = confirmed_df.trend_nonδ_lag00
    plot!(pl, ts_c, xs_nonδ_c; label="非デルタ株（推定・予測）",
        frange=zeros(Float64, nrow(confirmed_df)), fa=0.3, c=:cyan)
    plot!(pl, ts_c, xs_c; label="デルタ株（推定・予測）",
        frange=xs_nonδ_c, fa=0.3, c=:red3)
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_nonδ_c_l = confirmed_df[!,"trend_nonδ_"*lagname]
        if lag == 0
            plot!(pl, ts_c, xs_nonδ_c_l; label=:none, lc=:cyan)
        else
            plot!(pl, ts_c, xs_nonδ_c_l; label=:none, lc=:cyan, la=0.7)
        end        
    end
    plot!(pl, ts_c, xs_c; label=:none, lw=2, lc=:red3)
    plot!(pl, ts_c, cs_c; label=:none, m=:circle, lc=:royalblue, mc=:royalblue1)
    # 感染確認者予測
    ts_p = predict_df.date_value
    xs_p = predict_df.trend_lag00
    xs_nonδ_p = predict_df.trend_nonδ_lag00
    plot!(pl, ts_p, xs_nonδ_p; label=:none,
        frange=zeros(Float64, nrow(predict_df)), fa=0.1, c=:cyan, ls=:dot)
    plot!(pl, ts_p, xs_p; label=:none,
        frange=xs_nonδ_p, fa=0.1, c=:red3, lw=2, ls=:dot)
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_p_l = predict_df[!,"trend_"*lagname]
        xs_nonδ_p_l = predict_df[!,"trend_nonδ_"*lagname]
        if lag == 0
            plot!(pl, ts_p, xs_nonδ_p_l; label=:none, lc=:cyan, ls=:dot)
            plot!(pl, ts_p, xs_p_l; label=:none, lw=2, lc=:red3, ls=:dot)
        else
            plot!(pl, ts_p, xs_nonδ_p_l; label=:none, lc=:cyan, ls=:dot, la=0.7)
            plot!(pl, ts_p, xs_p_l; label=:none, lc=:red3, ls=:dot, la=0.7)
        end
    end
    # アノテーション
    annotate!(pl, date_to_value(plot_date_start) + 2, 100, text("""
        予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
        デルタ株確認者数の予測は最尤線をたどるとした場合""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left))
    plname = "confirmed"
    plots[plname] = pl
    @info plname
end
@info "--------"
@info "確認感染者数・トレンド（片対数グラフ）"
for once in 1:1
    pl = plot(
        size=(640,700),
#        size=(640,600),
        framestyle=:box,
        bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
        title="東京",
        legend=:topleft,
        fontfamily="Meiryo",
        yscale=:log
    )
    # x 軸設定
    sdate_value = date_to_value(plot_date_start)
    edate_value = date_to_value(plot_date_end)
    xlims!(pl, sdate_value, edate_value)
    nticks = length(date_ticks_dates)
    vdates = Float64[]
    labels = String[]
    for i in 1:nticks
        date = date_ticks_dates[i]
        if plot_date_start ≤ date ≤ plot_date_end
            push!(vdates, date_to_value(date))
            push!(labels, date_ticks_dict[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, 150.0, 8000.0)
#    ylims!(pl, 150.0, 4000.0)
    ylabel!(pl, "感染確認者数 [人/日]（対数目盛）")
    yticks!(pl,
        [150.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 4000.0, 6000.0, 8000.0],
        ["150", "200", "400", "600", "800", "1000", "1500", "2000", "4000", "6000", "8000"]
#        [150.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 4000.0],
#        ["150", "200", "400", "600", "800", "1000", "1500", "2000", "4000"]
    )
    # 感染確認者数
    ts_c = confirmed_df.date_value
    cs_c = confirmed_df.count
    xs_c = confirmed_df.trend
    plot!(pl, ts_c, cs_c; label="感染確認者数", m=:circle, lc=:royalblue, mc=:royalblue1)
    plot!(pl, ts_c, xs_c; label="トレンド", lw=2, lc=:red3)
    xs_nonδ_c = confirmed_df.trend_nonδ_lag00
    plot!(pl, ts_c, xs_nonδ_c; label="非デルタ株（推定・予測）",
        frange=zeros(Float64, nrow(confirmed_df)), fa=0.3, c=:cyan)
    plot!(pl, ts_c, xs_c; label="デルタ株（推定・予測）",
        frange=xs_nonδ_c, fa=0.3, c=:red3)
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_nonδ_c_l = confirmed_df[!,"trend_nonδ_"*lagname]
        if lag == 0
            plot!(pl, ts_c, xs_nonδ_c_l; label=:none, lc=:cyan)
        else
            plot!(pl, ts_c, xs_nonδ_c_l; label=:none, lc=:cyan, la=0.7)
        end        
    end
    plot!(pl, ts_c, xs_c; label=:none, lw=2, lc=:red3)
    plot!(pl, ts_c, cs_c; label=:none, m=:circle, lc=:royalblue, mc=:royalblue1)
    # 感染確認者予測
    ts_p = predict_df.date_value
    xs_p = predict_df.trend_lag00
    xs_nonδ_p = predict_df.trend_nonδ_lag00
    plot!(pl, ts_p, xs_nonδ_p; label=:none,
        frange=zeros(Float64, nrow(predict_df)), fa=0.1, c=:cyan, ls=:dot)
    plot!(pl, ts_p, xs_p; label=:none,
        frange=xs_nonδ_p, fa=0.1, c=:red3, lw=2, ls=:dot)
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_p_l = predict_df[!,"trend_"*lagname]
        xs_nonδ_p_l = predict_df[!,"trend_nonδ_"*lagname]
        if lag == 0
            plot!(pl, ts_p, xs_nonδ_p_l; label=:none, lc=:cyan, ls=:dot)
            plot!(pl, ts_p, xs_p_l; label=:none, lw=2, lc=:red3, ls=:dot)
            annotate!(pl, date_to_value(plot_date_end) - 2.0, xs_nonδ_p_l[end-14],
                text("変異株報告に\n追加の遅延がないとき", font("Meiryo",7), RGB(0.3,0.6,0.6), :right))
        else
            plot!(pl, ts_p, xs_nonδ_p_l; label=:none, lc=:cyan, ls=:dot, la=0.7)
            plot!(pl, ts_p, xs_p_l; label=:none, lc=:red3, ls=:dot, la=0.7)
            annotate!(pl, date_to_value(plot_date_end) - 2.0, xs_nonδ_p_l[end-14],
                text("追加の遅延 $lag 日", font("Meiryo",7), RGB(0.3,0.6,0.6), :right))
        end
    end
    # アノテーション
    annotate!(pl, date_to_value(plot_date_start) + 2, 170, text("""
        予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
        デルタ株確認者数の予測は最尤線をたどるとした場合""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left))
    plname = "confirmed_log"
    plots[plname] = pl
    @info plname
end
@info "--------"
@info "in English"
for once in 1:1
    pl = plot(
        size=(640,750),
        framestyle=:box,
        bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=32px,
        title="Tokyo",
        legend=:topleft,
        fontfamily="Meiryo",
        yscale=:log
    )
    # x 軸設定
    sdate_value = date_to_value(plot_date_start)
    edate_value = date_to_value(plot_date_end)
    xlims!(pl, sdate_value, edate_value)
    nticks = length(date_ticks_dates)
    vdates = Float64[]
    labels = String[]
    for i in 1:nticks
        date = date_ticks_dates[i]
        if plot_date_start ≤ date ≤ plot_date_end
            push!(vdates, date_to_value(date))
            push!(labels, date_ticks_dict_en[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, 150.0, 6000.0)
    ylabel!(pl, "Confirmed cases [cases/day] (logscale)")
    yticks!(pl,
        [150.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 4000.0, 6000.0],
        ["150", "200", "400", "600", "800", "1000", "1500", "2000", "4000", "6000"]
    )
    # 感染確認者数
    ts_c = confirmed_df.date_value
    cs_c = confirmed_df.count
    xs_c = confirmed_df.trend
    plot!(pl, ts_c, cs_c; label="Confirmed cases", m=:circle, lc=:royalblue, mc=:royalblue1)
    plot!(pl, ts_c, xs_c; label="Trend", lw=2, lc=:red3)
    xs_nonδ_c = confirmed_df.trend_nonδ_lag00
    plot!(pl, ts_c, xs_nonδ_c; label="non-δ cases (est./pred.)",
        frange=zeros(Float64, nrow(confirmed_df)), fa=0.3, c=:cyan)
    plot!(pl, ts_c, xs_c; label="δ cases (est./pred.)",
        frange=xs_nonδ_c, fa=0.3, c=:red3)
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_nonδ_c_l = confirmed_df[!,"trend_nonδ_"*lagname]
        if lag == 0
            plot!(pl, ts_c, xs_nonδ_c_l; label=:none, lc=:cyan)
        else
            plot!(pl, ts_c, xs_nonδ_c_l; label=:none, lc=:cyan, la=0.7)
        end        
    end
    plot!(pl, ts_c, xs_c; label=:none, lw=2, lc=:red3)
    plot!(pl, ts_c, cs_c; label=:none, m=:circle, lc=:royalblue, mc=:royalblue1)
    # 感染確認者予測
    ts_p = predict_df.date_value
    xs_p = predict_df.trend_lag00
    xs_nonδ_p = predict_df.trend_nonδ_lag00
    plot!(pl, ts_p, xs_nonδ_p; label=:none,
        frange=zeros(Float64, nrow(predict_df)), fa=0.1, c=:cyan, ls=:dot)
    plot!(pl, ts_p, xs_p; label=:none,
        frange=xs_nonδ_p, fa=0.1, c=:red3, lw=2, ls=:dot)
    for lag in likelihood_lags
        lagname = @sprintf("lag%02d", lag)
        xs_p_l = predict_df[!,"trend_"*lagname]
        xs_nonδ_p_l = predict_df[!,"trend_nonδ_"*lagname]
        if lag == 0
            plot!(pl, ts_p, xs_nonδ_p_l; label=:none, lc=:cyan, ls=:dot)
            plot!(pl, ts_p, xs_p_l; label=:none, lw=2, lc=:red3, ls=:dot)
            annotate!(pl, date_to_value(plot_date_end) - 2.0, xs_nonδ_p_l[end-14],
                text("no additional delay\nin reporting mutations", font("Meiryo",7), RGB(0.3,0.6,0.6), :right))
        else
            plot!(pl, ts_p, xs_nonδ_p_l; label=:none, lc=:cyan, ls=:dot, la=0.7)
            plot!(pl, ts_p, xs_p_l; label=:none, lc=:red3, ls=:dot, la=0.7)
            annotate!(pl, date_to_value(plot_date_end) - 2.0, xs_nonδ_p_l[end-14],
                text("$lag-day delay", font("Meiryo",7), RGB(0.3,0.6,0.6), :right))
        end
    end
    # アノテーション
    annotate!(pl, date_to_value(plot_date_start) + 2, 180, text("""
        The predictions are based on the expansion coefficient for the recent week
        and assume that the equivalent response would continue.
        The prediction of the δ cases assumes that the maximum likelihood line is followed.""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left))
    plname = "confirmed_log_en"
    plots[plname] = pl
    @info plname
end
@info "--------"
@info "拡大速度係数"
for lag in likelihood_lags
    lagname = @sprintf("lag%02d", lag)
    title = "東京"
    if lag != 0
        title *= "（変異株報告の追加の遅延 $(lag) 日）"
    end
    pl = plot(
        size=(640,600),
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
    for i in 1:nticks
        date = date_ticks_dates[i]
        if plot_date_start ≤ date ≤ plot_date_end
            push!(vdates, date_to_value(date))
            push!(labels, date_ticks_dict[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, -0.075, 0.150)
    ylabel!(pl, "拡大係数（時定数の逆数）[/日]")
    # グリッド
    plot!(pl, [sdate_value, edate_value], [0.0, 0.0]; label=:none, lc=RGB(0.3,0.3,0.3))
    # 拡大係数
    ts_c = confirmed_df.date_value
    ξds_c = confirmed_df.lgrowth
    ξds_nonδ_c = confirmed_df[!,"lgrowth_nonδ_"*lagname]
    plot!(pl, ts_c, ξds_c; label="拡大係数", lw=2, lc=:red3)
    plot!(pl, ts_c, ξds_nonδ_c; label="拡大係数（非デルタ株）推定", lc=:cyan)
    plot!(pl, ts_c, ξds_c; label=:none, lw=2, lc=:red3)
    # 拡大係数予測
    ts_p = predict_df.date_value
    ξds_p = predict_df[!,"lgrowth_"*lagname]
    ξds_nonδ_p = predict_df[!,"lgrowth_nonδ_"*lagname]
    plot!(pl, ts_p, ξds_nonδ_p; label="拡大係数（非デルタ株）想定値", lw=2, lc=:cyan, ls=:dot)
    plot!(pl, ts_p, ξds_p; label="デルタ株を加えた拡大係数予測", lw=2, lc=:red3, ls=:dot)
    # アノテーション
    annotate!(pl, date_to_value(plot_date_start) + 2, -0.06, text("""
        ここで拡大係数とは確認者数に指数関数をあてはめたときの時定数の逆数. 正で拡大、負で減少
        拡大係数が正のとき、その逆数の約 0.7 倍が拡大の倍加時間にあたる（0.05 なら倍加時間約 14 日）
        予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
        デルタ株拡大の予測は最尤線をたどるとした場合""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left))
    plname = "growth_"*lagname
    plots[plname] = pl
    @info plname
end
@info "--------"
@info "拡大速度係数（週あたり拡大率）"
for lag in likelihood_lags
    lagname = @sprintf("lag%02d", lag)
    title = "東京"
    if lag != 0
        title *= "（変異株報告の追加の遅延 $(lag) 日）"
    end
    pl = plot(
        size=(640,600),
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
    for i in 1:nticks
        date = date_ticks_dates[i]
        if plot_date_start ≤ date ≤ plot_date_end
            push!(vdates, date_to_value(date))
            push!(labels, date_ticks_dict[date])
        end
    end
    xticks!(pl, vdates, labels; xtickfontrotation=90.0)
    # y 軸設定
    ylims!(pl, -0.075, 0.150)
    ylabel!(pl, "週あたり拡大率")
    yt_wrs = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8]
    yticks!(pl,
        [log(wr)/7.0 for wr in yt_wrs],
        ["$(wr)" for wr in yt_wrs]
    )
    # グリッド
    plot!(pl, [sdate_value, edate_value], [0.0, 0.0]; label=:none, lc=RGB(0.3,0.3,0.3))
    # 拡大係数
    ts_c = confirmed_df.date_value
    ξds_c = confirmed_df.lgrowth
    ξds_nonδ_c = confirmed_df[!,"lgrowth_nonδ_"*lagname]
    plot!(pl, ts_c, ξds_c; label="拡大係数", lw=2, lc=:red3)
    plot!(pl, ts_c, ξds_nonδ_c; label="拡大係数（非デルタ株）推定", lc=:cyan)
    plot!(pl, ts_c, ξds_c; label=:none, lw=2, lc=:red3)
    # 拡大係数予測
    ts_p = predict_df.date_value
    ξds_p = predict_df[!,"lgrowth_"*lagname]
    ξds_nonδ_p = predict_df[!,"lgrowth_nonδ_"*lagname]
    plot!(pl, ts_p, ξds_nonδ_p; label="拡大係数（非デルタ株）想定値", lw=2, lc=:cyan, ls=:dot)
    plot!(pl, ts_p, ξds_p; label="デルタ株を加えた拡大係数予測", lw=2, lc=:red3, ls=:dot)
    # アノテーション
    annotate!(pl, date_to_value(plot_date_start) + 2, -0.06, text("""
        ここで拡大係数とは確認者数に指数関数をあてはめたときの時定数の逆数. 正で拡大、負で減少
        拡大係数が正のとき、その逆数の約 0.7 倍が拡大の倍加時間にあたる（0.05 なら倍加時間約 14 日）
        予測は、直近 1 週間の拡大係数をもとに同等の対応が継続する（非デルタ株に対し一定）と想定
        デルタ株拡大の予測は最尤線をたどるとした場合""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left))
    plname = "growth_m_"*lagname
    plots[plname] = pl
    @info plname
end

#=
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
=#

end #module
