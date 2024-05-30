module output_save_functions

using XLSX, CSV, DataFrames, StatsBase, Dates

include("basic_functions.jl")
include("choice_file.jl")



#--- Save aggregated reults to Excel sheets
function save_aggregated_results_as_excel(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    
    OUT_TIMESERIES = preprocess_aggregated_output(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    OUT_OPTIMIZATION = preprocess_aggregated_OPTINFO(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)


    # check if results folder exists, else create it
    str_folder = abspath(joinpath(pwd(), "../results"))
    if !(isdir(str_folder))
        mkdir(str_folder)
    end

    # Save OVERVIEW on results in Excel sheet
    str_file = string(str_folder *  "/" * Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * " - OVERVIEW on Output Rolling Optimization.xlsx")
    XLSX.writetable(str_file, "Timeseries" => OUT_TIMESERIES, "Optimization" => OUT_OPTIMIZATION)
    print("Saved results successfully - file: " * str_folder * "\n")


end

#--- Save multiple matrices and vectors of the reults to Excel-file
function save_results(results, results_names)
    # Define path and file name as string
    str_path = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * " - Output Rolling Optimization.xlsx"
    # Create excel-sheet
    XLSX.openxlsx(str_path, mode="w") do xf
        for i in 1:length(results)
            if i == 1
                sheet = xf[1]
                XLSX.rename!(sheet, results_names[i])
            else
                sheet = XLSX.addsheet!(xf, results_names[i])
            end
            sheet["A1"] = convert(Array, transpose(results[i]))
        end
    end

    print("Saved results successfully - file: " * str_path * "\n")
end

function quantiles_results(MAT_RES, MAT_WEIGHT, TIME, WEIGHT)
    # Initialize output
    eval_res = zeros(size(MAT_RES,2),5)
    # Get minimum values
    eval_res[:,1] = minimum(MAT_RES, dims = 1)'
    # Get 25 % quantile
    eval_res[:,2] = [quantile(MAT_RES[:,t], 0.25) for t in TIME]
    if WEIGHT
        # Get weighted average values
        eval_res[:,3] = sum(MAT_RES .* (MAT_WEIGHT) , dims = 1)' ./ sum(MAT_WEIGHT , dims = 1)'
    else
        # Get average values
        eval_res[:,3] = [quantile(MAT_RES[:,t], 0.5) for t in TIME]
    end
    # Get 75 % quantile
    eval_res[:,4] = [quantile(MAT_RES[:,t], 0.75) for t in TIME]
    # Get maximum values
    eval_res[:,5] = maximum(MAT_RES, dims = 1)'

    return eval_res
end

#--- Get quantiles of power flows accross AC lines
function quantiles_line_loading(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    # Get index of timeseries input as temporary calc time
    TMP_CALC_TIME = collect(1:length(PARAMETER_SETTINGS.CALC_TIME))
    # Number of quantiles
    N = 3
    # Get set of AC lines
    SET_AC_LINES = MODEL_INPUT.SETS.BRANCHES[size(INPUT.IN_TRAFOS,1)+1:end]
    # Get absolute values of relative power flows, i.e. normalized using the rated power
    ABS_FLOW_EL_LINE = broadcast(abs, RESULTS.FLOW_EL_BRANCH[SET_AC_LINES,:] ./ MODEL_INPUT.BRANCHES.PMAX[SET_AC_LINES])
    # Initialize matrix for results of quantile analysis
    QUANT_ABS_FLOW_EL_LINE = zeros(length(TMP_CALC_TIME),N+2)
    # Get names of quantiles as string matrix for plot
    quant_names = reshape([string(trunc(Int,q/(N+1)*100)) * " %" for q in 1:N],1,N)
    # Get 1 % and 99 % quantiles
    QUANT_ABS_FLOW_EL_LINE[:,1] = [quantile(ABS_FLOW_EL_LINE[:,t], 0.01) for t in TMP_CALC_TIME ]
    QUANT_ABS_FLOW_EL_LINE[:,N+2] = [quantile(ABS_FLOW_EL_LINE[:,t], 0.99) for t in TMP_CALC_TIME ]
    # For loop, to get quantile values
    for q in 1:N
        QUANT_ABS_FLOW_EL_LINE[:,q+1] = [quantile(ABS_FLOW_EL_LINE[:,t], q/(N+1)) for t in TMP_CALC_TIME ]
    end
    return QUANT_ABS_FLOW_EL_LINE
end

#--- Pre-process output
function preprocess_aggregated_output(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
#-- Convert datatype of prices from Matrix with Missing values to Matrix of Float64
    # Get index of timeseries input as temporary calc time
    TMP_CALC_TIME = collect(1:length(PARAMETER_SETTINGS.CALC_TIME))
    # Get prices as matrix
    PRICES = Matrix{Float64}(RESULTS.PRICES)
    # Get weights as matrix
    MAT_WEIGHT = Matrix(MODEL_INPUT.ASSIGNMENTS.A_L) * MODEL_INPUT.LOADS.TIMESERIES'
    # Get quantiles of prices
    eval_prices = quantiles_results(PRICES, MAT_WEIGHT, TMP_CALC_TIME, true)

    # Get quantiles of brnach loading
    eval_load_flow = quantiles_line_loading(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)

    # Get absolute values of HVDC power flows
    ABS_FLOW_HVDC = broadcast(abs, RESULTS.HVDC_FLOW)
    # Get HVDC loading
    HVDC_LOADING = Vector{Float64}(Matrix{Float64}(sum(ABS_FLOW_HVDC, dims = 1)')[:,1])
    # Get relative HVDC loading
    HVDC_REL_LOADING = Vector{Float64}(Matrix{Float64}(sum(ABS_FLOW_HVDC ./ MODEL_INPUT.HVDCS.P_MAX, dims = 1)')[:,1]) ./ MODEL_INPUT.SETS.HVDCS[end]
    replace!(HVDC_REL_LOADING, NaN=>0)

#-- Aggregate results
    # Get aggregated filling level of storage
    TMP_FILLING_LEVEL = sum(RESULTS.FILLING_LEVEL, dims = 1)'
    # Get aggregated net power injections of storages due to charging
    AGG_CHARGE_STOR = Vector{Float64}(Matrix{Float64}(sum(RESULTS.STORAGES_CHARGING, dims = 1)')[:,1])
    # Get aggregated net power injections of storages due to discharging
    AGG_DISCHARGE_STOR = Vector{Float64}(Matrix{Float64}(sum(RESULTS.STORAGES_DISCHARGING, dims = 1)')[:,1])
    # Get aggregated net power injections of storages
    AGG_GEN_STOR = AGG_DISCHARGE_STOR - AGG_CHARGE_STOR
        # Vector{Float64}([TMP_FILLING_LEVEL[idx+1]-TMP_FILLING_LEVEL[idx] for idx in 1:size(TMP_FILLING_LEVEL,1)-1])
    # Get aggregated net power injections of PV
    AGG_GEN_PV = Vector{Float64}(sum(MODEL_INPUT.RENEWABLES.TIMESERIES[:,1:size(INPUT.IN_PV,1)], dims = 2)[:,1])
    # Get aggregated net power injections of Wind
    AGG_GEN_WIND = Vector{Float64}(sum(MODEL_INPUT.RENEWABLES.TIMESERIES[:,size(INPUT.IN_PV,1)+1:end], dims = 2)[:,1])
    # Get aggregated net power injections of hydro run of the river
    AGG_GEN_HY_ROR = Vector{Float64}(sum(MODEL_INPUT.HY_ROR.TIMESERIES, dims = 2)[:,1])
    # Get curtailed EE generation
    AGG_RENEWABLE_CURTAILED = Vector{Float64}(sum(RESULTS.RENEWABLE_FEED_IN', dims = 2)[:,1])
    # Get aggregated renewable curtailment
    AGG_RENEWABLE_CURTAILMENT = AGG_GEN_WIND + AGG_GEN_PV + AGG_GEN_HY_ROR - AGG_RENEWABLE_CURTAILED
    # Get aggregated net power injections of electricity only units
    AGG_GEN_EL_ONLY = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_EL_ONLY, dims = 1)')[:,1])
    # Get aggregated net power injections of CHP units
    AGG_GEN_EL_CHP = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_EL_CHP, dims = 1)')[:,1])
    # Get aggregated net power injections of BIO units
    AGG_GEN_EL_BIO = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_EL_BIO, dims = 1)')[:,1])
    # Get aggregated net power demand of loads
    AGG_LOAD = Vector{Float64}(sum(MODEL_INPUT.LOADS.TIMESERIES, dims = 2)[:,1])
    # Heat pumps
    AGG_HEATPUMP = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_HEAT_PUMP, dims = 1)')[:,1])
    # Electric mobility
    AGG_LOAD_EMOB_RESIDENTIAL = Vector{Float64}(sum(RESULTS.LOAD_EMOB_RESIDENTIAL', dims = 2)[:,1])
    AGG_LOAD_EMOB_COMMERCIAL = Vector{Float64}(sum(RESULTS.LOAD_EMOB_COMMERCIAL', dims = 2)[:,1])
    AGG_LOAD_EMOB_QUICKCHARGE = Vector{Float64}(sum(RESULTS.LOAD_EMOB_QUICKCHARGE', dims = 2)[:,1])


    # Get aggregated net power injections of loads
    AGG_SLACKS_POS = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_EL_SLACKS_POS, dims = 1)')[:,1])
    AGG_SLACKS_NEG = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_EL_SLACKS_NEG, dims = 1)')[:,1])

    # Get aggregated net power injections of CHP units
    AGG_GEN_HEAT_CHP = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_HEAT_CHP, dims = 1)')[:,1])
    # Get aggregated heat production of power to heat units
    AGG_GEN_HEAT_PTH = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_HEAT_PTH, dims = 1)')[:,1])
    # Get aggregated heat production of heat only units
    AGG_GEN_HEAT_ONLY = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_HEAT_ONLY, dims = 1)')[:,1])
    # Get total aggregated net heat demand
    AGG_HEAT_DEMAND_ALL = Vector{Float64}(sum(MODEL_INPUT.HEAT_DEMAND.TIMESERIES, dims = 2)[:,1]) .+ sum(MODEL_INPUT.CONST_HEAT_LOAD.Q_MAX)
    # Get timeseries-based aggregated net heat demand
    AGG_HEAT_DEMAND_TIMESERIES = Vector{Float64}(sum(MODEL_INPUT.HEAT_DEMAND.TIMESERIES, dims = 2)[:,1])
    # Get aggregated const. heat load
    AGG_HEAT_DEMAND_CONST_HEAT_LOAD = AGG_HEAT_DEMAND_ALL - AGG_HEAT_DEMAND_TIMESERIES
    # Get heat generation of slacks
    AGG_SLACK_Q = Vector{Float64}(Matrix{Float64}(sum(RESULTS.SLACK_Q, dims = 1)')[:,1])
    # Get aggregated hydrogen production of power to gas units
    AGG_GEN_HYDROGEN = Vector{Float64}(RESULTS.GEN_PTG' * MODEL_INPUT.PTG.ETA_P_MAX)
    # Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_HYDROGEN, dims = 1)')[:,1])
    # Get aggregated hydrogen production of power to gas units
    AGG_GEN_PTG = Vector{Float64}(Matrix{Float64}(sum(RESULTS.GEN_PTG, dims = 1)')[:,1])

    # Get total aggregated net power demand of loads and heatpumps
    AGG_LOAD_ALL = AGG_LOAD + AGG_HEATPUMP + AGG_GEN_HEAT_PTH + AGG_GEN_PTG + AGG_LOAD_EMOB_RESIDENTIAL + AGG_LOAD_EMOB_COMMERCIAL + AGG_LOAD_EMOB_QUICKCHARGE

#-- Consolidate and save results in excel sheets
    # Consolidate all aggregated results in one Matrix
    AGG_RES = mapreduce(permutedims, vcat,
        [AGG_LOAD_ALL,
        AGG_LOAD,
        AGG_HEATPUMP,
        AGG_GEN_HEAT_PTH,
        AGG_GEN_PTG,
        AGG_LOAD_EMOB_RESIDENTIAL,
        AGG_LOAD_EMOB_COMMERCIAL,
        AGG_LOAD_EMOB_QUICKCHARGE,
        AGG_GEN_PV,
        AGG_GEN_WIND,
        AGG_GEN_HY_ROR,
        AGG_RENEWABLE_CURTAILED,
        AGG_RENEWABLE_CURTAILMENT,
        AGG_GEN_EL_BIO,
        AGG_GEN_EL_ONLY,
        AGG_GEN_EL_CHP,
        AGG_CHARGE_STOR,
        AGG_DISCHARGE_STOR,
        AGG_GEN_STOR,
        HVDC_LOADING,
        HVDC_REL_LOADING,
        AGG_SLACKS_POS,
        AGG_SLACKS_NEG,
        AGG_HEAT_DEMAND_ALL,
        AGG_HEAT_DEMAND_TIMESERIES,
        AGG_HEAT_DEMAND_CONST_HEAT_LOAD,
        AGG_GEN_HEAT_CHP,
        AGG_GEN_HEAT_PTH,
        AGG_GEN_HEAT_ONLY,
        AGG_SLACK_Q,
        AGG_GEN_HYDROGEN])'
    # Convert Matrix with consolidated results to DataFrame
    OUT_AGG_RES = DataFrame([PARAMETER_SETTINGS.CALC_TIME AGG_RES eval_prices eval_load_flow], :auto)

    # Rename columns of DataFrame
    col_names = ["Timestep",
        "Load (el. all)",
        "Load (conv.)",
        "Load (heat pumps)",
        "Load (PtH)",
        "Load (PtG)",
        "Load (Emob residential)",
        "Load (Emob commercial)",
        "Load (Emob quickcharge)",
        "PV",
        "Wind",
        "Hydro (run-of-river)",
        "RES infeed after curtailment", # previous name "Curtailed renewables"
        "RES curtailed", # previous name: Curtailment
        "Biomass/Biogas",
        "Gas power plants",
        "CHP (el. gen.)",
        "Storages charging",
        "Storages discharging",
        "Storages net power injection",
        "HVDC abs.",
        "HVDC abs. (in %)",
        "Slacks_pos",
        "Slacks_neg",
        "Heat demand (all)",
        "Heat demand (district heating)",
        "Heat demand (industrial)",
        "CHP (th. gen.)",
        "Power-to-Heat",
        "Heat plant",
        "Heat slacks",
        "Power-to-Gas - Hydrogen prod.",
        "Min. LMPs",
        "25%-Quantile of LMPs",
        "AVG LMPs",
        "75%-Quantile of LMPs",
        "Max. LMPs",
        "Min. line loading",
        "25% quantile line loading",
        "AVG line loading",
        "75% quantile line loading",
        "Max. line loading"]
    rename!(OUT_AGG_RES, Symbol.(col_names))

    return OUT_AGG_RES
end


function preprocess_aggregated_OPTINFO(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)

    OUT_OPTIMIZIATION = DataFrame

    MODELSTATUS = Vector{String}(RESULTS.MODELSTATUS)
    OBJ_VAL = Vector{Float64}(RESULTS.OBJ_VAL)
    SOLVE_TIME = Vector{Float64}(RESULTS.SOLVE_TIME)
    
    SLACK_EL_POS = Vector{Float64}(RESULTS.AGG_SLACK_EL_POS)
    SLACK_EL_NEG = Vector{Float64}(RESULTS.AGG_SLACK_EL_NEG)
    SLACK_HP = Vector{Float64}(RESULTS.AGG_SLACK_HP)
    SLACK_Q = Vector{Float64}(RESULTS.AGG_SLACK_Q)
 


    OUT_OPTIMIZIATION = hcat(MODELSTATUS, OBJ_VAL, SOLVE_TIME, SLACK_EL_POS, SLACK_EL_NEG, SLACK_HP, SLACK_Q)
    DF_SHEET = DataFrame(OUT_OPTIMIZIATION, :auto)
    rename!(DF_SHEET,[:SOVLESTATUS,:OBJ_VAL,:SOLVE_TIME, :SLACK_EL_POS, :SLACK_EL_NEG, :SLACK_HP, :SLACK_Q]) # rename column names

    return DF_SHEET
end

end