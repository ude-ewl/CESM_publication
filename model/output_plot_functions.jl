module output_plot_functions

using Plots
using Statistics
using ArchGDAL
using CSV
using DataFrames
using CategoricalArrays
#using Gadfly

function plot_selection(input_dict, MODEL_INPUT, RESULTS_ALL, RD_RESULTS_ALL, PARAMETER_SETTINGS)

    if input_dict["LMP_aggregated"]
        if PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL+RD"
            plot_lmp_aggregated(MODEL_INPUT, RD_RESULTS_ALL)
        else
            plot_lmp_aggregated(MODEL_INPUT, RESULTS_ALL)
        end
    end

    if input_dict["nodal_prices"]
        plot_nodal_prices(input_dict["date"], PARAMETER_SETTINGS)
    end

    if input_dict["line_loading"]
        plot_line_loading(input_dict["date"], PARAMETER_SETTINGS)
    end

    if input_dict["rd_quantities"]
        plot_rd_quantities(input_dict["date"], input_dict["rd_date"], PARAMETER_SETTINGS)
    end

    return nothing
end


function plot_lmp_aggregated(MODEL_INPUT, RD_RESULTS_ALL)

    cols_prices = ["Max. LMP", "95% Percentile", "75%-Quantile", "Median LMP", "25%-Quantile", "5% Percentile", "Min. LMP"]
    
    size(RD_RESULTS_ALL.PRICES)
    data_sorted = sort(RD_RESULTS_ALL.PRICES, dims=1)
    relevant_rows_zero = round.(Int, size(RD_RESULTS_ALL.PRICES, 1)*[1, 0.95, 0.75, 0.5, 0.25, 0.05, 0.0])
    relevant_rows = relevant_rows_zero
    relevant_rows[end] = 1

    data_vals = data_sorted[relevant_rows,:]'

    bla = plot(data_vals, linetype=:steppre, labels=permutedims(cols_prices),legend=:best, title="LMP", ylabel="€/MWh", xlabel="hour")
    display(bla)

end

function plot_nuts1()
    filepath = "Y:\\Materialien\\Daten\\GEo\\NUTS3\\NUTS_2021\\NUTS_RG_01M_2021_4326_LEVL_1.shp"
    dataset = ArchGDAL.read(filepath)
    layer = ArchGDAL.getlayer(dataset, 0)
    numfeature = ArchGDAL.nfeature(layer)
    p = Plots.plot(size=(1200,1200), aspect_ratio=1.4, label=nothing, grid=nothing)
    for ct_region in 0:(numfeature-1)
        country_code, fs = ArchGDAL.getfeature(layer, ct_region) do f
                    return ArchGDAL.getfield(f, 2), ArchGDAL.getgeom(f)
                end
        if country_code == "DE"
            plot!(p,fs, color=:white)
        end
    end
    return p
end

function plot_nodal_prices(date, PARAMETER_SETTINGS)

    p = plot_nuts1()

    # ---- model input data ----
    lmps = DataFrame(CSV.File( abspath("../results/" * date * " - Results//EL_LMPs.csv")))
    lmps = lmps[:,"Timestep" .!= names(lmps)]
    lmps_trans = Matrix(lmps)'

    # take out "Line_-" from line name and take mean value as utilization value
    lmps_aggr = DataFrame( name=[x[6:end] for x in names(lmps)], util=mean.(eachcol(lmps)))

    # calculate new aggragated rows
    insertcols!(lmps_aggr, :mean_val => vec(mean(lmps_trans, dims=2)))

    # add hours
    for ct_hr in collect(1:convert(Int,length(PARAMETER_SETTINGS.CALC_TIME)/PARAMETER_SETTINGS.PERIOD_LENGTH))
        colname = "h" * string(convert(Int,ct_hr))
        lmps_aggr[!, colname] = Array(lmps[ct_hr,:])
    end
   
    elec_nodes = DataFrame(CSV.File(abspath("../tmp_data/elec_network_nodes.csv")))
    elec_locations = DataFrame(CSV.File(abspath("../tmp_data/location.csv")))
    elec_nodes = leftjoin(elec_nodes, elec_locations[:,["location_id", "location_name", "latitude", "longitude"]], on="location_id")
    print(string(100*sum(ismissing.(elec_nodes.latitude))/length(elec_nodes.latitude)) * " % missing latlons due to location_id mismatch")

    # plot LMPS
    prices = lmps_aggr[!,"mean_val"]
    scatter!(elec_nodes.longitude, elec_nodes.latitude, marker_z=prices, 
                xaxis=false, yaxis=false, legend = false, right_margin = 10Plots.mm,
                colorbar = true, color=:jet, colorbar_title = " \nLMP in €/MWh")
    display(p)

end


function plot_line_loading(date, PARAMETER_SETTINGS)

    highVoltageThreshold = 110
    highload_threshold = 0.95

    # ---- static model input data ----
    elec_nodes = DataFrame(CSV.File(abspath("../tmp_data/elec_network_nodes.csv")))
    elec_ac_lines = DataFrame(CSV.File(abspath("../tmp_data/elec_network_ac_lines.csv")))
    elec_ac_lines = elec_ac_lines[elec_ac_lines.voltage_level.>highVoltageThreshold,:] # only hv lines
    elec_dc_lines = DataFrame(CSV.File(abspath("../tmp_data/elec_network_hvdc_lines.csv")))
    elec_locations = DataFrame(CSV.File(abspath("../tmp_data/location.csv")))
    # join locations into elec_nodes
    elec_nodes = leftjoin(elec_nodes, elec_locations[:,["location_id", "location_name", "latitude", "longitude"]], on="location_id")

    # ---- line loading results ----
    lineload_ac = DataFrame(CSV.File( abspath("../results/" * date * " - Results/EL_AC Lines (rel.).csv")))
    lineload_ac = lineload_ac[:,"Timestep" .!= names(lineload_ac)]
    lineload_ac_trans = Matrix(lineload_ac)'
    lineload_dc = DataFrame(CSV.File( abspath("../results/" * date * " - Results/EL_HVDC Flow.csv")))


    # ---- obtain/aggregate hourly loadings to single value that should be plotted
    # ---- select one row to use ----
    #values_to_plot = mean.(eachcol(lineload))

    highload = convert(Matrix{Float64},Matrix(lineload_ac.>highload_threshold))
    values_to_plot = vec(sum(highload, dims=1)')

    # take out "Line_-" from line name and take mean value as utilization value
    lineload_ac_aggr = nothing
    lineload_ac_aggr = DataFrame( name=[x[6:end] for x in names(lineload_ac)], util=values_to_plot)

    lineload_ac_aggr[ismissing.(lineload_ac_aggr[!,"util"]), "util"] .= NaN   # replace missings by nan, since missings cannot be plotted
    lineload_aggr = lineload_ac_aggr
    # join line loading into lines dataframe
    elec_ac_lines = leftjoin(elec_ac_lines, lineload_aggr, on="name")
    elec_ac_lines[ismissing.(elec_ac_lines[!,"util"]), "util"] .= NaN

    # ---- check consistency of data before geoplot ----
    #p_hist = plot()
    #histogram(lineload_aggr.util)
    #maximum(lineload_aggr.util)
    #plot(lineload_aggr.util)

    p = plot_nuts1() # plot background shapes


    function plotLineLoading(elec_nodes, elec_lines, col_name, color_seperators, colorlist, legend_multiplier, legend_prefix, legend_suffuix)

        # for debugging
        lines_w_values = 0
        lines_wo_values = 0
        vec_lon = []
        vec_lat = []
        vec_line = []
        vec_util = []
    
        col = elec_lines[:,col_name]
    
        for ct in eachrow(elec_lines)
    
            # detect start and end coordinate
            idx1 = collect(elec_nodes[!,"node_id"] .== ct.node_id1)
            idx2 = collect(elec_nodes[!,"node_id"] .== ct.node_id2)
    
            # if found
            if idx1 != nothing && idx2 != nothing
    
                lat1 =  elec_nodes[ collect(idx1),["latitude"]].latitude[1]
                lon1 =  elec_nodes[ collect(idx1),["longitude"]].longitude[1]
                lat2 =  elec_nodes[ collect(idx2),["latitude"]].latitude[1]
                lon2 =  elec_nodes[ collect(idx2),["longitude"]].longitude[1]
    
                # if there are coordinates and plotting is thus possible, stor information
                vec_lon = append!(vec_lon, [NaN, lon1, lon2])
                vec_lat = append!(vec_lat, [NaN, lat1, lat2])
                vec_line = append!(vec_line, [NaN, ct.ac_line_id, ct.ac_line_id ])
                vec_util = append!(vec_util, [NaN, ct[col_name], ct[col_name]])
    
                # for debugging
                lines_w_values = lines_w_values + 1
            else
                lines_wo_values = lines_wo_values + 1
            end
        end
    
        # debug plot for line loading
        #plot(vec_util[.!(isnan.(vec_util))])
    
        # all lines in gray
        plot!(vec_lon, vec_lat, color=:lightgray, labels="all", 
            xaxis=false, yaxis=false, )
    
        threshold_increasing = color_seperators[1] < color_seperators[end]
    
        # highlight lines with high load and draw above
        ct_color = 1
        for ct_threshold in color_seperators
            vec_util_mod = vec_util
            vec_util_mod[isnan.(vec_util_mod)] .= 9999999
            id_show = vec_util_mod.>ct_threshold
            vec_lon = vec_lon[id_show]
            vec_lat = vec_lat[id_show]
            vec_line = vec_line[id_show]
            vec_util = vec_util[id_show]
            plot!(vec_lon, vec_lat, width=3, color=colorlist[ct_color], labels=string(legend_prefix, string(convert(Int,legend_multiplier*ct_threshold)), legend_suffix))
            ct_color += 1
        end
        plot!(aspect_ratio=1, grid=0,  xlims=(5,19.5)) # xaxis=0, yaxis=0, xtick=0, ytick=0,
    
    end

    color_seperators = [50, 100, 150, 200, 250]
    colorlist        = ["gold", "orange", "orangered", "red", "darkred"]
    legend_multiplier = 1  # 100 if percent
    legend_prefix = "> " * string(convert(Int64,100*highload_threshold)) * "% for "
    legend_suffix = " h"
    plotLineLoading(elec_nodes, elec_ac_lines, "util", color_seperators, colorlist, legend_multiplier, legend_prefix, legend_suffix)
    display(p)
    #savefig("lineloading_aggregated.svg")

end




function plot_rd_quantities(date_market, date_rd, PARAMETER_SETTINGS)

    p = plot_nuts1()

    technology_files =  Dict( 
        "EL_Gen CHP.csv"=>+1,
        "EL_Gen Elec Only.csv" =>+1,
        "EL_Hydro-ROR.csv"=>0,
        "EL_Load (HP).csv"=>-1,
        "EL_Load (PTG).csv"=>-1,
        "EL_Load (PTH).csv"=>-1,
        "EL_Pumped-storage charging.csv TECH:BS"=>-1,
        "EL_Pumped-storage charging.csv TECH:PS"=>-1,
        "EL_Pumped-storage discharging.csv TECH:BS"=>+1,
        "EL_Pumped-storage discharging.csv TECH:PS"=>+1,
        "EL_Renewable curtailed (PV, Wind & Hydro-ROR).csv"=>+1,
        "EL_EMOB_residential.csv"=>+1,
        "EL_EMOB_comercial.csv"=>+1,
        "EL_EMOB_quickcharge.csv"=>+1)

    schedule_df = DataFrame("name" => [])

    for (filestring, sign) in technology_files
        filestring_short = split(filestring, " TECH:")[1]
        if sign != 0
            market_df = DataFrame(CSV.File( abspath("../results/" * date_market * " - Results//" * filestring_short)))  # gibt es immer in beiden ergebnissen die gleichen Spalten?!?!?
            rd_df = DataFrame(CSV.File( abspath("../results/" * date_rd * " - Results//" * filestring_short)))

            market = market_df[:,"Timestep" .!= names(market_df)] # check if same columns
            rd = rd_df[:,"Timestep" .!= names(rd_df)]

            if (sum(names(market) .== names(rd)) < length(names(market)))
                println("columns in market result and rd result do not match for " * filestring)
                errorABC
            end

            # filter storage types
            if (filestring == "EL_Pumped-storage charging.csv TECH:BS") | (filestring == "EL_Pumped-storage discharging.csv TECH:PS")
                market_df[:,occursin.("_BS_", names(market_df))]
                rd_df[:,occursin.("_BS_", names(rd_df))]
            elseif (filestring == "EL_Pumped-storage ischarging.csv TECH:BS") | (filestring == "EL_Pumped-storage discharging.csv TECH:PS")
                market_df[:,occursin.("_PS_", names(market_df))]
                rd_df[:,occursin.("_PS_", names(rd_df))]
            end
            
            # get difference between results and aggregate over all time steps
            rd_quantities = sign * (Matrix(rd) - Matrix(market))
            df_temp = DataFrame(rd_quantities, names(market))
            mean_val = mean.(eachcol(df_temp))

            if mean(mean_val) == 0
                println(string(filestring * " no schedule change, technology is skipped"))
            else
                # transpose information
                rd_quantities_df_transposed = DataFrame([[names(df_temp)]; collect.(eachrow(df_temp))], [:column; Symbol.(axes(df_temp, 1))]) # transpose
                headerstring = replace(string(filestring), ".csv" => "")
                rd_quantities_df_transposed[:,headerstring] = mean_val # add additional column with mean value
                rd_quantities_df_transposed.name = [s[1:12] for s in names(market)] # add additional column with mean value
                df_short = rd_quantities_df_transposed[:,["name",headerstring]]
                df_nonneg = df_short[ abs.(df_short[:,2]) .>= 0.001,:]
                # sum up values at same elec_location
                numcols = names(df_nonneg, findall(x -> eltype(x) <: Number, eachcol(df_nonneg)))
                df_unique = combine(groupby(df_nonneg, ["name"]), numcols .=> sum .=> numcols)

                # merge technology rd in technology-matrix
                schedule_df = outerjoin(schedule_df, df_unique, on=:name)
            end
        end
    end
    schedule_df = coalesce.(schedule_df, 0.0) # fill missing with zero

    println(names(schedule_df))

    # add all technologies
    schedule_df[!,"nodal_sum"] = vec(sum(Matrix(schedule_df[:,2:end]),dims=2))
    _, temp_tuples = findmax(Matrix(schedule_df[:,2:end-1]), dims=2)
    col_info = getindex.(temp_tuples[:,1],2)
    
    dict_techs = Dict(1 => "EL_Gen CHP", 2 => "EL_Gen Elec Only", 3 => "EL_Load (PTG)", 4 => "EL_Load (PTH)", 5 => "EL_Renewable curtailed (PV, Wind & Hydro-ROR)", 6 => "EL_Load (HP)", 7 => "EL_Pumped-storage charging BS", 8 => "EL_Pumped-storage charging PS", 9 => "EL_Pumped-storage discharging BS", 10 => "EL_Pumped-storage discharging PS")
    dict_color = Dict(1 => "red",        2 => "sienna",           3 => "purple",        4 => "gray", 5 => "green", 6 => "gold2", 7 => "cyan", 8 => "blue", 9 => "cyan", 10 => "blue")

    rd_techs = [dict_color[col_value] for col_value in col_info]
    #println(rd_techs)
   
    #schedule_df[!,"major_rd_tech"] = categorical(rd_techs)
    schedule_df[!,"major_rd_tech"] = rd_techs


    #println(schedule_df[!,"major_rd_tech"])
    #println(schedule_df)

    # load node to location assginement 
    elec_nodes = DataFrame(CSV.File(abspath("../tmp_data/elec_network_nodes.csv")))
    names_short = [s[1:12] for s in schedule_df.name]
    nodes_idx = in.(elec_nodes.name , [names_short])

    # load coordiantes from locations
    elec_locations = DataFrame(CSV.File(abspath("../tmp_data/location.csv")))
    elec_nodes = leftjoin(elec_nodes, elec_locations[:,["location_id", "location_name", "latitude", "longitude"]], on="location_id")
 
    println(string(100*sum(ismissing.(elec_nodes.latitude))/length(elec_nodes.latitude)) * " % missing latlons due to location_id mismatch")

    elec_nodes = leftjoin(elec_nodes, schedule_df, on="name")
    
    idx_nonmissing = .!(ismissing.(elec_nodes.nodal_sum))


    scatter!(elec_nodes[idx_nonmissing, "longitude"], elec_nodes[idx_nonmissing, "latitude"], 
                alpha=0.5, xaxis=false, yaxis=false, label = nothing, right_margin = 10Plots.mm,
                color = elec_nodes[idx_nonmissing, "major_rd_tech"], title = " \n Nodal aggregated redispatch quantities for optimization period",
                markershape = :utriangle, markersize=elec_nodes[idx_nonmissing, "nodal_sum"]/5, markerstrokewidth=0)

    # annotate max value            
    elec_nodes_nonzero = coalesce.(elec_nodes, 0.0) # fill missing with zero
    annotate!(14,49.5,string(round(maximum(elec_nodes_nonzero[:,:nodal_sum]),digits=1)) * " MWh", :left)
    max_id = findall(elec_nodes_nonzero[:,:nodal_sum].==maximum(elec_nodes_nonzero[:,:nodal_sum]))
    plot!(LinRange(elec_nodes[max_id,"longitude"][1],13.5,100),LinRange(elec_nodes[max_id,"latitude"][1],49.3,100), color="black", label=nothing) # line to largest circle
    plot!(LinRange(13.5,16,100),LinRange(49.3,49.3,100), color="black", label=nothing, xaxis=false, yaxis=false) # horiz line

    # annotate min value            
    annotate!(14,50.2,string(round(minimum(elec_nodes_nonzero[:,:nodal_sum]),digits=1)) * " MWh", :left)
    min_id = findall(elec_nodes_nonzero[:,:nodal_sum].==minimum(elec_nodes_nonzero[:,:nodal_sum]))
    plot!(LinRange(elec_nodes[min_id,"longitude"][1],13.5,100),LinRange(elec_nodes[min_id,"latitude"][1],50,100), color="black", label=nothing) # line to largest circle
    plot!(LinRange(13.5,16,100),LinRange(50,50,100), color="black", label=nothing, xaxis=false, yaxis=false) # horiz line

    display(p)

    println("Offset: " * string(sum(schedule_df.nodal_sum)) * " MW")

end

# ==== unused stuff =====


#plot([sum(MODEL_INPUT.RENEWABLES.TIMESERIES, dims = 2),
#sum(RESULTS_ALL.GEN_EL_ONLY, dims = 1)',
#sum(RESULTS_ALL.GEN_EL_CHP, dims = 1)',
#sum(RESULTS_ALL.GEN_EL_BIO, dims = 1)',
# sum(RESULTS.GEN_EL_HY_ROR, dims = 1)',
#sum(MODEL_INPUT.LOADS.TIMESERIES, dims = 2) +
#sum(MODEL_INPUT.HEATPUMPS.ELEC_TIMESERIES, dims = 2)],
#linetype=:steppre, labels=["RE" "EL_Only" "CHP" "BIO" "HY_ROR" "Load & Heat pump demand"], legend=:right)


# # cols_prices = ["Timestep", "Min. LMPs", "25%-Quantile of LMPs", "AVG LMPs", "75%-Quantile of LMPs", "Max. LMPs"]
# prices = Matrix{Float64}(OUT_AGG_RES[!,cols_prices])
#
# plot(Matrix{Float64}(prices[!,cols_prices])[:,1], Matrix{Float64}(OUT_AGG_RES[!,cols_prices])[:,2],
#     fillrange = 1, fillalpha = 0.35, c = 1, label = "Confidence band")
#
# sum(sum(RD_RESULTS_ALL.GEN_EL_ONLY,dims=1)' - sum(RESULTS_ALL.GEN_EL_ONLY,dims=1)')
# sum(sum(RD_RESULTS_ALL.FILLING_LEVEL,dims=1)' - sum(RESULTS_ALL.FILLING_LEVEL,dims=1)')
# sum(sum(RD_RESULTS_ALL.PRICES,dims=1)' - sum(RESULTS_ALL.PRICES,dims=1)')
# sum(sum(RD_RESULTS_ALL.GEN_EL_SLACKS_POS,dims=1)' - sum(RESULTS_ALL.GEN_EL_SLACKS_POS,dims=1)')
#
# bla1 = sum(RD_RESULTS_ALL.RENEWABLE_FEED_IN,dims=1)'
# bla2 = sum(RESULTS_ALL.RENEWABLE_FEED_IN,dims=1)'
# plot(bla1)
#
#
# market_rd = DataFrame(variable=[], type=[], market=[], rd=[], market_minus_rd=[])
#
# # sum
# vars_sum = [:OPERATIONAL_COSTS, :GEN_EL_ONLY, :GEN_EL_CHP, :GEN_HEAT_CHP, :GEN_EL_BIO, :RENEWABLE_FEED_IN, :GEN_PTG, :GEN_HYDROGEN, :GEN_HEAT_PTH, :GEN_HEAT_ONLY, :GEN_HEAT_PUMP]
# for ct_var in vars_sum
#     row_data = DataFrame(variable = ct_var, type="SUM", market = sum(RESULTS_ALL[ct_var]), rd = sum(RD_RESULTS_ALL[ct_var]), market_minus_rd = sum(RESULTS_ALL[ct_var])-sum(RD_RESULTS_ALL[ct_var]))
#     append!(market_rd, row_data)
# end
#
# # mean
# vars_mean = [:PRICES, :FILLING_LEVEL, :HVDC_FLOW, :FLOW_EL_BRANCH]
# for ct_var in vars_mean
#     row_data = DataFrame(variable = ct_var, type="MEAN", market = mean(RESULTS_ALL[ct_var]), rd = mean(RD_RESULTS_ALL[ct_var]), market_minus_rd = mean(RESULTS_ALL[ct_var])-mean(RD_RESULTS_ALL[ct_var]))
#     append!(market_rd, row_data)
# end
#
# # max
# vars_max = [:GEN_EL_ONLY, :GEN_EL_CHP, :GEN_EL_BIO, :RENEWABLE_FEED_IN, :GEN_PTG, :GEN_HYDROGEN, :GEN_HEAT_PTH, :GEN_HEAT_PUMP]
# for ct_var in vars_max
#     row_data = DataFrame(variable = ct_var, type="MAX", market = mean(RESULTS_ALL[ct_var]), rd = mean(RD_RESULTS_ALL[ct_var]), market_minus_rd = mean(RESULTS_ALL[ct_var])-mean(RD_RESULTS_ALL[ct_var]))
#     append!(market_rd, row_data)
# end
#
#
# pwd()
# CSV.write("market_rd_comparison.csv", market_rd, delim=";", decimal=',')
#
# # gen elec only
# plot( Array[ sum(RESULTS_ALL.GEN_EL_ONLY,dims=1)', sum(RD_RESULTS_ALL.GEN_EL_ONLY,dims=1)'], label=["market" "rd"], title="GEN_EL_ONLY" )
#
# # renewable infeed
# plot( Array[ sum(RESULTS_ALL.RENEWABLE_FEED_IN,dims=1)', sum(RD_RESULTS_ALL.RENEWABLE_FEED_IN,dims=1)'], label=["market" "rd"], title="RENEWABLE_FEED_IN" )


#---
end

