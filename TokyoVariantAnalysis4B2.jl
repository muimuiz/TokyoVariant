module TokyoVariantAnalysis4B2
using Logging
using DataFrames
using Dates

using ..TokyoVariantAnalysis4: date_to_value, value_to_date
using ..TokyoVariantAnalysis4A

export heatmap_array
export hm_x_min, hm_x_max, hm_y_min, hm_y_max, hm_x_size, hm_y_size, hm_f

@info "module", @__MODULE__

@info "--------"
@info "ヒートマップ形式尤度設定"
const hm_x_size = 500
const hm_y_size = 500
const hm_x_min = date_to_value(Date("2021-05-25"))
const hm_x_max = hm_x_min + 100.0
const hm_y_min = 0.0
const hm_y_max = 1.0
const hm_z_na = TokyoVariantAnalysis4A.lh_logL_threshold
@info "(hm_x_size, hm_y_size)", (hm_x_size, hm_y_size)
@info "(hm_x_min, hm_x_max)", (hm_x_min, hm_x_max)
@info "(hm_y_min, hm_y_max)", (hm_y_min, hm_y_max)
@info "hm_z_na", hm_z_na
heatmap_array = fill(hm_z_na, (hm_x_size, hm_y_size))
for once in 1:1
    b_array = falses(hm_x_size, hm_y_size)
    x_range = range(hm_x_min, hm_x_max; length=hm_x_size)
    for (i, row) in enumerate(eachrow(likelihoods_df))
        if i % 1000 == 0
            print(stderr, "\r$(round(Int, 100.0*i/nrow(likelihoods_df)))%")
        end
        iy_p::Union{Int,Nothing} = nothing
        for (ix, x) in enumerate(x_range)
            y = 1.0 / (1.0 + 1.0 / exp(row.α + row.β * x))
            iy = round(Int, (hm_y_size - 1.0) * (y - hm_y_min) / (hm_y_max - hm_y_min)) + 1
            if 1 ≤ iy ≤ hm_y_size
                rlogL = row.logL - likelihood_max.logL
                if !isnothing(iy_p) && (1 ≤ iy_p ≤ hm_y_size) && (iy_p != iy)
                    for iy_c in iy_p:sign(iy-iy_p):iy
                        if !b_array[ix, iy_c]
                            heatmap_array[ix, iy_c] = rlogL 
                            b_array[ix, iy_c] = true
                        end
                    end
                else
                    if !b_array[ix, iy]
                        heatmap_array[ix, iy] = rlogL
                        b_array[ix, iy] = true
                    end
                end
            end
            iy_p = iy
        end
    end
    print(stderr, "\r100%\n")
end
#
function hm_f(x, y)
    ix = round(Int, (hm_x_size - 1.0) * (x - hm_x_min) / (hm_x_max - hm_x_min)) + 1
    iy = round(Int, (hm_y_size - 1.0) * (y - hm_y_min) / (hm_y_max - hm_y_min)) + 1
    if (1 ≤ ix ≤ hm_x_size) && (1 ≤ iy ≤ hm_y_size)
        return exp(heatmap_array[ix, iy])
    else
        return exp(hm_z_na)
    end
end

end #module
