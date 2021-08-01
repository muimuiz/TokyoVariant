module TokyoVariantAnalysis4C2
using Logging
using DataFrames
using Dates
using Plots
using Plots.PlotMeasures

using ..TokyoVariantAnalysis4
using ..TokyoVariantAnalysis4A
using ..TokyoVariantAnalysis4B2

@info "--------"
@info "プロット・バックエンド・設定"
pyplot(
    titlefont=font("Meiryo",10),
    guidefont=font("Meiryo",10),
    tickfont=font("Meiryo",9)
)
const plot_date_start = Date("2021-05-25")
const plot_date_end   = plot_date_start + Day(100)
const date_ticks_dict = Dict(
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
@info "L452R 割合尤度プロット"
for once in 1:1
    pl = plot(
        size=(640,600),
        framestyle=:box,
        bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=-32px,
        title="東京（健安研調査・民間検査機関等調査）"
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
    ylims!(pl, 0.0, 1.0)
    ylabel!(pl, "L452R 変異のある変異株（主に δ 株）の割合")
    yticks!(pl,
        [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        ["0\\%", "20\\%", "40\\%", "60\\%", "80\\%", "100\\%"]
    )
    # 尤度ヒートマップ描画
    hm_xrange = range(hm_x_min, hm_x_max; length=hm_x_size)
    hm_yrange = range(hm_y_min, hm_y_max; length=hm_y_size)    
    heatmap!(pl, hm_xrange, hm_yrange, hm_f;
        c=:Wistia, colorbar_title="相対尤度（最大尤度を 1 とする）", colorbar_titlefontfamily="Meiryo")
    plot!(pl, (t)->1.0/(1.0+1.0/exp(likelihood_max.α+likelihood_max.β*t)), lc=:red3, ls=:dot)
    # データ点
    pname = "L452R_prop"
    lmtl = -100.0
    lmtu = +100.0
    lmt(r) = min(max(r, lmtl), lmtu)
    ds  = variants_df.date_value
    #des = df.date_error
    rs  = variants_df[!,Symbol("$(pname)")]
    ls  = variants_df[!,Symbol("$(pname)_cil")]
    us  = variants_df[!,Symbol("$(pname)_ciu")]
    n_rows = nrow(variants_df)
    pm_ds  = Vector{Float64}(); sizehint!(pm_ds,  n_rows)
    #pm_des = Vector{Float64}(); sizehint!(pm_des, n_rows)
    pm_rs  = Vector{Float64}(); sizehint!(pm_rs,  n_rows)
    pm_els = Vector{Float64}(); sizehint!(pm_els, n_rows)
    pm_eus = Vector{Float64}(); sizehint!(pm_eus, n_rows)
    ps_ds  = Vector{Float64}(); sizehint!(ps_ds,  n_rows)
    #ps_des = Vector{Float64}(); sizehint!(ps_des, n_rows)
    ps_rs  = Vector{Float64}(); sizehint!(ps_rs,  n_rows)
    ps_els = Vector{Float64}(); sizehint!(ps_els, n_rows)
    ps_eus = Vector{Float64}(); sizehint!(ps_eus, n_rows)
    for i in 1:n_rows
        if !ismissing(rs[i])
            r  = lmt(rs[i])
            el = r - lmt(ls[i])
            eu = lmt(us[i]) - r
            if TokyoVariantAnalysis4A.lh_sdate_value ≤ ds[i] ≤ TokyoVariantAnalysis4A.lh_edate_value
                push!(pm_ds,  ds[i])
                #push!(pm_des, des[i])
                push!(pm_rs,  r)
                push!(pm_els, el)
                push!(pm_eus, eu)
            else
                push!(ps_ds,  ds[i])
                #push!(ps_des, des[i])
                push!(ps_rs,  r)
                push!(ps_els, el)
                push!(ps_eus, eu)
            end
        end
    end
    scatter!(pl, pm_ds, pm_rs; yerrors=(pm_els, pm_eus),
        color=RGB(0.0,0.5,1.0), msc=RGB(0.4,0.4,0.4), legend=:none)
    scatter!(pl, ps_ds, ps_rs; yerrors=(ps_els, ps_eus),
        color=RGB(0.0,0.5,1.0), msc=RGB(0.7,0.7,0.7), legend=:none)
    # グリッド
    vline!(pl, vdates; lc=RGBA(0.0,0.0,0.0,0.1), ls=:dot)
    hline!(pl, [0.2, 0.4, 0.6, 0.8]; lc=RGBA(0.0,0.0,0.0,0.1), ls=:dot)
    # アノテーション
    #coef = ratio_fit_coefs[pname][2]
    annotate!(pl, date_to_value(Date("2021-05-28")), 0.925, text("""
        「東京都新型コロナウイルス感染症対策本部報」より
        民間検査機関等における L452R の検査が報告されるようになった 2021-06-01 以降のデータ
        ただし、デルタ株クラスターが含まれているとみられる場合は、クラスターあたり 1 人分のみ計上
        赤点線は指数関数モデルでの最尤線. 背景の色は最大尤度に相対的な尤度を表す""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left)
    )
    plname = "L452R"
    plots[plname] = pl
    @info plname
end
@info "--------"
@info "in English"
for once in 1:1
    pl = plot(
        size=(640,600),
        framestyle=:box,
        bottom_margin=8px, left_margin=16px, top_margin=16px, right_margin=-32px,
        title="Tokyo (Inst. of Public Health and Private Labs.)"
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
    ylims!(pl, 0.0, 1.0)
    ylabel!(pl, "Proportion of the L452R mutation (mainly the δ variant)")
    yticks!(pl,
        [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        ["0\\%", "20\\%", "40\\%", "60\\%", "80\\%", "100\\%"]
    )
    # 尤度ヒートマップ描画
    hm_xrange = range(hm_x_min, hm_x_max; length=hm_x_size)
    hm_yrange = range(hm_y_min, hm_y_max; length=hm_y_size)    
    heatmap!(pl, hm_xrange, hm_yrange, hm_f;
        c=:Wistia, colorbar_title="relative likelihood (the maximum likelihood is 1)", colorbar_titlefontfamily="Meiryo")
    plot!(pl, (t)->1.0/(1.0+1.0/exp(likelihood_max.α+likelihood_max.β*t)), lc=:red3, ls=:dot)
    # データ点
    pname = "L452R_prop"
    lmtl = -100.0
    lmtu = +100.0
    lmt(r) = min(max(r, lmtl), lmtu)
    ds  = variants_df.date_value
    #des = df.date_error
    rs  = variants_df[!,Symbol("$(pname)")]
    ls  = variants_df[!,Symbol("$(pname)_cil")]
    us  = variants_df[!,Symbol("$(pname)_ciu")]
    n_rows = nrow(variants_df)
    pm_ds  = Vector{Float64}(); sizehint!(pm_ds,  n_rows)
    #pm_des = Vector{Float64}(); sizehint!(pm_des, n_rows)
    pm_rs  = Vector{Float64}(); sizehint!(pm_rs,  n_rows)
    pm_els = Vector{Float64}(); sizehint!(pm_els, n_rows)
    pm_eus = Vector{Float64}(); sizehint!(pm_eus, n_rows)
    ps_ds  = Vector{Float64}(); sizehint!(ps_ds,  n_rows)
    #ps_des = Vector{Float64}(); sizehint!(ps_des, n_rows)
    ps_rs  = Vector{Float64}(); sizehint!(ps_rs,  n_rows)
    ps_els = Vector{Float64}(); sizehint!(ps_els, n_rows)
    ps_eus = Vector{Float64}(); sizehint!(ps_eus, n_rows)
    for i in 1:n_rows
        if !ismissing(rs[i])
            r  = lmt(rs[i])
            el = r - lmt(ls[i])
            eu = lmt(us[i]) - r
            if TokyoVariantAnalysis4A.lh_sdate_value ≤ ds[i] ≤ TokyoVariantAnalysis4A.lh_edate_value
                push!(pm_ds,  ds[i])
                #push!(pm_des, des[i])
                push!(pm_rs,  r)
                push!(pm_els, el)
                push!(pm_eus, eu)
            else
                push!(ps_ds,  ds[i])
                #push!(ps_des, des[i])
                push!(ps_rs,  r)
                push!(ps_els, el)
                push!(ps_eus, eu)
            end
        end
    end
    scatter!(pl, pm_ds, pm_rs; yerrors=(pm_els, pm_eus),
        color=RGB(0.0,0.5,1.0), msc=RGB(0.4,0.4,0.4), legend=:none)
    scatter!(pl, ps_ds, ps_rs; yerrors=(ps_els, ps_eus),
        color=RGB(0.0,0.5,1.0), msc=RGB(0.7,0.7,0.7), legend=:none)
    # グリッド
    vline!(pl, vdates; lc=RGBA(0.0,0.0,0.0,0.1), ls=:dot)
    hline!(pl, [0.2, 0.4, 0.6, 0.8]; lc=RGBA(0.0,0.0,0.0,0.1), ls=:dot)
    # アノテーション
    #coef = ratio_fit_coefs[pname][2]
    annotate!(pl, date_to_value(Date("2021-05-28")), 0.925, text("""
        Data source: Daily reports of “Tōkyō-to Shingata Corona Virus Kansenshō Taisaku Hombu”
        Using data since 2021-06-01 when private laboratories began reporting L452R tests.
        Only one case included per cluster if case numbers contain Delta clusters.
        The red dotted line is the maximum likelihood line for the exponential model.
        The background color indicates the likelihood relative to the maximum.""",
        font("Meiryo",7), RGB(0.3,0.3,0.3), :left)
    )
    plname = "L452R_en"
    plots[plname] = pl
    @info plname
end

end # module
