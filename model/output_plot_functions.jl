module output_plot_functions

using Plots
using Statistics
using ArchGDAL
using CSV
using DataFrames
using CategoricalArrays
#using Gadfly

function plot_selection(input_dict, MODEL_INPUT, RESULTS_ALL, PARAMETER_SETTINGS)

    if input_dict["LMP_aggregated"]
        plot_lmp_aggregated(MODEL_INPUT, RESULTS_ALL)
        
    end

    if input_dict["nodal_prices"]
        plot_nodal_prices(input_dict["date"], PARAMETER_SETTINGS)
    end

    if input_dict["line_loading"]
        plot_line_loading(input_dict["date"], PARAMETER_SETTINGS)
    end

    return nothing
end


function plot_lmp_aggregated(MODEL_INPUT, RESULTS_ALL)

    cols_prices = ["Max. LMP", "95% Percentile", "75%-Quantile", "Median LMP", "25%-Quantile", "5% Percentile", "Min. LMP"]
    
    size(RESULTS_ALL.PRICES)
    data_sorted = sort(RESULTS_ALL.PRICES, dims=1)
    relevant_rows_zero = round.(Int, size(RESULTS_ALL.PRICES, 1)*[1, 0.95, 0.75, 0.5, 0.25, 0.05, 0.0])
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
   
    elec_nodes = DataFrame(CSV.File(abspath("../input_data/elec_network_nodes.csv")))
    elec_locations = DataFrame(CSV.File(abspath("../input_data/location.csv")))
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
    elec_nodes = DataFrame(CSV.File(abspath("../input_data/elec_network_nodes.csv")))
    elec_ac_lines = DataFrame(CSV.File(abspath("../input_data/elec_network_ac_lines.csv")))
    elec_ac_lines = elec_ac_lines[elec_ac_lines.voltage_level.>highVoltageThreshold,:] # only hv lines
    elec_dc_lines = DataFrame(CSV.File(abspath("../input_data/elec_network_hvdc_lines.csv")))
    elec_locations = DataFrame(CSV.File(abspath("../input_data/location.csv")))
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

end

