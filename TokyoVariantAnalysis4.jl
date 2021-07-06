module TokyoVariantAnalysis4
using Logging
using CSV
using DataFrames
using Dates
using HypothesisTests: BinomialTest, confint

export variants_df
export date_to_value, value_to_date

@info "module", @__MODULE__

function has_name(df::DataFrame, name::String)
    return !isnothing(findfirst(isequal(name), names(df)))
end

@info "--------"
@info "作業ディレクトリ"
@info pwd()

@info "--------"
@info "基準日"
const epoch_date = Date("2021-01-01")
@info "epoch_date", epoch_date
#
date_to_value(d::Date) = Dates.value(d - epoch_date)
value_to_date(v)::Date = epoch_date + Day(v)

@info "========"
@info "データ読み込み・変換"

const organization_names = ["k", "p"]   # k 健安研, p 民間検査機関等
const variant_names = ["N501Y", "L452R"]
@info "organization_names", organization_names
@info "variant_names", variant_names

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
    @info "日付 (start_date, end_date, period_days, mid_date, date_value, date_error)"
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
    @info "--------"
    @info "変異株検査総数のコピー (ORG_n)"
    for orgname in organization_names
        colname = orgname*"_n"  # カラム名
        @assert has_name(csv_df, colname)
        df[!,colname] = csv_df[!,colname]
    end
    @info size(df)
    @info "--------"
    @info "変異株カウントのコピー (ORG_VAR)"
    for orgname in organization_names
        for varname in variant_names
            colname = orgname*"_"*varname   # カラム名
            @assert has_name(csv_df, colname)
            df[!,colname] = csv_df[!,colname]
        end
    end
    @info size(df)

    n_rows = nrow(df)
    @info "--------"
    @info "検査機関にわたる変異株検査総数の合計 (n)"
    ns::Vector{Union{Int,Missing}} = fill(missing, n_rows)
    for i in 1:n_rows
        n = missing
        for orgname in organization_names
            colname = orgname*"_n"  # 検査機関別カラム名
            orgn = df[i,colname]
            if !ismissing(orgn)
                if ismissing(n)
                    n  = orgn
                else
                    n += orgn
                end
            end
        end
        ns[i] = n
    end
    df.n = ns
    @info size(df)
    @info "--------"
    @info "検査機関にわたる変異株カウント合計 (VAR)"
    for varname in variant_names
        cs::Vector{Union{Int,Missing}} = fill(missing, n_rows)
        for i in 1:n_rows
            c = missing
            for orgname in organization_names
                colname = orgname*"_"*varname   # 検査機関別カラム名
                orgc = df[i,colname]
                if !ismissing(orgc)
                    if ismissing(c)
                        c  = orgc
                    else
                        c += orgc
                    end
                end
            end
            cs[i] = c
        end
        df[!,varname] = cs
    end
    @info size(df)

    @info "--------"
    @info "クラスター・渡航関連者を補正した L452R 変異株数と検査総数"
    @assert has_name(csv_df, "L452R_tr")   # 渡航者接触者
    @assert has_name(csv_df, "L452R_cn")   # クラスター数
    @assert has_name(csv_df, "L452R_cp")   # クラスター内人数
    ns_adj = copy(df.n)
    cs_adj = copy(df.L452R)
    for i in 1:n_rows
        n_adj = ns_adj[i]
        c_adj = cs_adj[i]
        tr = csv_df.L452R_tr[i]
        if !ismissing(tr)
            n_adj -= tr
            c_adj -= tr
        end
        cn = csv_df.L452R_cn[i]
        cp = csv_df.L452R_cp[i]
        if !ismissing(cp)
            n_adj -= cp
            c_adj -= cp
            if !ismissing(cn)
                n_adj += cn
                c_adj += cn
            end
        end
        @assert ismissing(n_adj) || (n_adj ≥ 0)
        @assert ismissing(c_adj) || (c_adj ≥ 0)
        ns_adj[i] = n_adj
        cs_adj[i] = c_adj
    end
    df.n     = ns_adj
    df.L452R = cs_adj
    return df
end

const tokyo_variants_csv_filename = "東京変異株4.csv"
variants_df = load_df(tokyo_variants_csv_filename)

@info "========"
@info "統計量算出"

const variant_ratio_pairs = [
    ("L452R", "N501Y"),
]

function statistics!(df::DataFrame)
    n_rows = nrow(df)
    @info "--------"
    @info "全体に対する割合を計算 (VAR_prop)"
    for varname in variant_names
        df[!,varname*"_prop"] = df[!,varname] ./ df.n
    end
    @info size(df)
    @info "--------"
    @info "1 変異株の全体に対する信頼区間を計算 (VAR_prop_cil, VAR_prop_ciu)"
    for varname in variant_names
        cs = df[!,varname]
        ns = df.n
        cils::Vector{Union{Float64,Missing}} = fill(missing, n_rows)
        cius::Vector{Union{Float64,Missing}} = fill(missing, n_rows)
        for i in 1:n_rows
            (cils[i], cius[i]) = if !ismissing(cs[i]) && !ismissing(ns[i])
                confint(BinomialTest(cs[i], ns[i]))
            else
                (missing, missing)
            end
        end
        df[!,varname*"_prop_cil"] = cils
        df[!,varname*"_prop_ciu"] = cius
    end
    @info size(df)
    @info "--------"
    @info "1 変異株のオッズを計算 (VAR_odds)"
    for varname in variant_names
        df[!,varname*"_odds"] = 1.0 ./ ((df.n ./ df[!,varname]) .- 1.0)
    end
    @info size(df)
    @info "--------"
    @info "オッズの信頼区間を計算 (VAR_odds_cil, VAR_odds_ciu)"
    for varname in variant_names
        cs = df[!,varname]
        ns = df.n
        cils::Vector{Union{Float64,Missing}} = fill(missing, n_rows)
        cius::Vector{Union{Float64,Missing}} = fill(missing, n_rows)
        for i in 1:n_rows
            (cils[i], cius[i]) = if !ismissing(cs[i]) && !ismissing(ns[i])
                1.0 ./ (1.0 ./ confint(BinomialTest(cs[i], ns[i])) .- 1.0)
            else
                (missing, missing)
            end
        end
        df[!,varname*"_odds_cil"] = cils
        df[!,varname*"_odds_ciu"] = cius
    end
    @info size(df)
    @info "--------"
    @info "オッズの対数（ロジット） (log_VAR_odds, log_VAR_odds_cil, log_VAR_odds_ciu)"
    for varname in variant_names
        oddname = varname*"_odds"
        df[!,"log_"*oddname]        = log.(df[!,oddname])
        df[!,"log_"*oddname*"_cil"] = log.(df[!,oddname*"_cil"])
        df[!,"log_"*oddname*"_ciu"] = log.(df[!,oddname*"_ciu"])
    end
    @info "--------"
    @info "特定 2 変異株間の比を計算 (VAR1_VAR2)"
    for (varname1, varname2) in variant_ratio_pairs
        @info varname1, varname2
        df[!,varname1*"_"*varname2] = df[!,varname1] ./ df[!,varname2]
    end
    @info size(df)
    @info "--------"
    @info "2 変異株間の比の信頼区間を計算 (VAR1_VAR2_cil, VAR1_VAR2_ciu)"
    @info "2 変異株のみからなる 2 項分布による近似"
    for (varname1, varname2) in variant_ratio_pairs
        @info varname1, varname2
        cs1 = df[!,varname1]
        cs2 = df[!,varname2]
        cils::Vector{Union{Float64,Missing}} = fill(missing, n_rows)
        cius::Vector{Union{Float64,Missing}} = fill(missing, n_rows)
        for i in 1:n_rows
            (cils[i], cius[i]) = if !ismissing(cs1[i]) && !ismissing(cs2[i])
                1.0 ./ (1.0 ./ confint(BinomialTest(cs1[i], cs1[i] + cs2[i])) .- 1.0)
            else
                (missing, missing)
            end
        end
        df[!,varname1*"_"*varname2*"_cil"] = cils
        df[!,varname1*"_"*varname2*"_ciu"] = cius
    end
    @info size(df)
    @info "--------"
    @info "2 変異株間比の対数 (log_VAR1_VAR2, log_VAR1_VAR2_cil, log_VAR1_VAR2_ciu)"
    for (varname1, varname2) in variant_ratio_pairs
        @info varname1, varname2
        pairname = varname1*"_"*varname2
        df[!,"log_"*pairname]        = log.(df[!,pairname])
        df[!,"log_"*pairname*"_cil"] = log.(df[!,pairname*"_cil"])
        df[!,"log_"*pairname*"_ciu"] = log.(df[!,pairname*"_ciu"])
    end
    @info size(df)
end

statistics!(variants_df)

@info "--------"
@info "解析データ書き出し"
const output_csv_filename = "tokyo_variant_analysis4.csv"
@info output_csv_filename
CSV.write(output_csv_filename, variants_df; missingstring="NA")
@info size(variants_df)

end #module
