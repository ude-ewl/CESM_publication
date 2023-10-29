module output_save_functions

using Dates, XLSX, DataFrames, CSV, StatsBase #, JLD2

include("basic_functions.jl")
include("choice_file.jl")


function save_output_combination(RESULTS_ALL, MODEL_INPUT, INPUT, COLUMN_DEFINITIONS, input_dict, PARAMETER_SETTINGS)

    # Save aggregated results
    if input_dict["aggregated_xlsx"]
        save_aggregated_results_as_excel(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    end

    # Save detailed results as CSVs
    if input_dict["detailed_csv"]
        save_detailed_results_as_csv(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS, COLUMN_DEFINITIONS)
    end

    #Save detailed results for online model
    if input_dict["online_model"]
        save_detailed_results_for_online_model(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS, COLUMN_DEFINITIONS)
    end

    # Save AoTRF results
    if input_dict["AoTRF"]
        save_results_for_AoTRF(MODEL_INPUT,RESULTS_ALL)
    end

    # Save as JLD2
    if input_dict["jld2"]
        filenames = ["PARAMETER_SETTINGS", "INPUT", "MODEL_INPUT", "RESULTS_ALL"];
        save_multiple_jld2([PARAMETER_SETTINGS, INPUT, MODEL_INPUT, RESULTS_ALL],  filenames)
          
    end

    return nothing
end


#--- Save reults to Excel-file
function save_results_as_excel_ffs(RESULTS)
    # Define path and file name as string
    str_folder = abspath(joinpath(pwd(), "../results"))
    str_path = str_folder * "/" * Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * " - Output Rolling Optimization.xlsx"
    XLSX.openxlsx(str_path, mode="w") do xf
        for i in 1:length(RESULTS)
            if i == 1
                sheet = xf[1]
                XLSX.rename!(sheet, String.(collect(keys(RESULTS)))[i])
            else
                sheet = XLSX.addsheet!(xf, String.(collect(keys(RESULTS)))[i])
            end
            sheet["A1"] = convert(Array, transpose(RESULTS[i]))
        end
    end
end

function save_results_as_excel_ffs_2(RESULTS, NAME)
    # Define path and file name as string
    str_folder = abspath(joinpath(pwd(), "../results"))
    str_path = str_folder * "/" * Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * "-" * NAME *".xlsx"
    XLSX.openxlsx(str_path, mode="w") do xf
        for i in 1:length(RESULTS)
            if i == 1
                sheet = xf[1]
                XLSX.rename!(sheet, String.(collect(keys(RESULTS)))[i])
            else
                sheet = XLSX.addsheet!(xf, String.(collect(keys(RESULTS)))[i])
            end
            if typeof(RESULTS[i]) != Matrix{Union{Missing, Float64}}
                sheet["A1"] = convert(Array, RESULTS[i])
            else
                sheet["A1"] = convert(Array, transpose(RESULTS[i]))
            end
        end
    end
end

#--- Function to detailed results as CSV files
function save_detailed(IN_DF, col_names, time, sub_folder, file_name)
    # Get column names
    tmp_cols = broadcast(string,vcat("Timestep",col_names))
    # Get data as matrix
    tmp_mat = hcat(time, Matrix{Float64}(IN_DF))
    # Create DataFrame
    tmp_df = DataFrame(tmp_mat, tmp_cols)
    # Save data to csv
    CSV.write(sub_folder*"\\"*file_name*".csv", tmp_df)
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

#--- Save reults to Excel-file
function save_aggregated_results_as_excel(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    
    OUT_AGG_RES = preprocess_aggregated_output(RESULTS_ALL, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    
    # check if results folder exists, else create it
    str_folder = abspath(joinpath(pwd(), "../results"))
    if !(isdir(str_folder))
        mkdir(str_folder)
    end
    # Save OVERVIEW on results in Excel
    str_file = string(str_folder *  "/" * Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * " - OVERVIEW on Output Rolling Optimization.xlsx")
    XLSX.writetable(str_file, OVERVIEW=( collect(DataFrames.eachcol(OUT_AGG_RES)), DataFrames.names(OUT_AGG_RES) ))
    print("Saved AGGREGATED results successfully - file: " * str_folder * "\n")
end

# Save detailed results as csv files
function save_detailed_results_as_csv(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS, COLUMN_DEFINITIONS)
    # Get index of timeseries input as temporary calc time
    TMP_CALC_TIME = collect(1:length(PARAMETER_SETTINGS.CALC_TIME))
    # Define folder name
    folder_path = abspath(joinpath(pwd(), "../results/")) * Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * " - Results"
    # Create folder if not exists
    if !(isdir(folder_path))
        mkdir(folder_path)
    end
    sleep(10)
    # Save electrical load
    save_detailed(MODEL_INPUT.LOADS.TIMESERIES, names(INPUT.IN_LOAD_TIMESERIES)[2:end],
        TMP_CALC_TIME, folder_path, "EL_Load (conv.)")
    # Save heat pump's load
    save_detailed(RESULTS.GEN_HEAT_PUMP', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, folder_path, "EL_Load (HP)")
    # Save heat pump's load
    save_detailed(RESULTS.GEN_HEAT_PTH', INPUT.IN_PTH_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_PTH.unit_name],
        TMP_CALC_TIME, folder_path, "EL_Load (PTH)")
    # Save heat pump's load
    save_detailed(RESULTS.GEN_PTG', INPUT.IN_PTG[:,COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.unit_name],
        TMP_CALC_TIME, folder_path, "EL_Load (PTG)")
    # Save PV generation
    save_detailed(MODEL_INPUT.RENEWABLES.TIMESERIES[:,1:size(INPUT.IN_PV,1)], names(INPUT.IN_RENEWABLE_TIMESERIES)[2:size(INPUT.IN_PV,1)+1],
        TMP_CALC_TIME, folder_path, "EL_PV")
    # Save Wind generation
    save_detailed(MODEL_INPUT.RENEWABLES.TIMESERIES[:,size(INPUT.IN_PV,1)+1:end], names(INPUT.IN_RENEWABLE_TIMESERIES)[2+size(INPUT.IN_PV,1):end],
        TMP_CALC_TIME, folder_path, "EL_Wind")
    # Save hydro ron-of-the-river generation
    save_detailed(MODEL_INPUT.HY_ROR.TIMESERIES, names(INPUT.IN_HY_ROR_TIMESERIES)[2:end],
        TMP_CALC_TIME, folder_path, "EL_Hydro-ROR")
    # Save curtailed renewable feed-in
    save_detailed(RESULTS.RENEWABLE_FEED_IN', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, folder_path, "EL_Renewable curtailed (PV, Wind & Hydro-ROR)")
    # Save biomass/biogas generation
    save_detailed(RESULTS.GEN_EL_BIO', INPUT.IN_BIO[:,COLUMN_DEFINITIONS.COLUMNS_BIO.unit_name],
        TMP_CALC_TIME, folder_path, "EL_Biomass & Biogas")
    # Save electricty only generation
    save_detailed(RESULTS.GEN_EL_ONLY', INPUT.IN_CONVENTIONAL[:,COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.unit_name],
        TMP_CALC_TIME, folder_path, "EL_Gen Elec Only")
    # Save electricty only generation costs
    save_detailed(RESULTS.OPERATIONAL_COSTS', INPUT.IN_CONVENTIONAL[:,COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.unit_name],
        TMP_CALC_TIME, folder_path, "Costs Elec Only")
    # Save CHP generation
    save_detailed(RESULTS.GEN_EL_CHP', INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.unit_name],
        TMP_CALC_TIME, folder_path, "EL_Gen CHP")
    # Save CHP generation costs
    save_detailed(RESULTS.OPERATIONAL_COSTS_CHP', INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.unit_name],
    TMP_CALC_TIME, folder_path, "Costs CHP")
    # Save pumped storages charging
    save_detailed(RESULTS.STORAGES_CHARGING', INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.unit_name],
        TMP_CALC_TIME, folder_path, "EL_Storage charging")
    # Save pumped storages discharging
    save_detailed(RESULTS.STORAGES_DISCHARGING', INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.unit_name],
        TMP_CALC_TIME, folder_path, "EL_Storage discharging")
    # Save HVDC flow
    save_detailed(RESULTS.HVDC_FLOW', INPUT.IN_HVDC[:,COLUMN_DEFINITIONS.COLUMNS_HVDC.unit_name],
        TMP_CALC_TIME, folder_path, "EL_HVDC Flow")
    # Save pumped storages discharging
    save_detailed(RESULTS.GEN_EL_SLACKS_NEG', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, folder_path, "EL_Slacks Generation")
    # Save pumped storages charging
    save_detailed(RESULTS.GEN_EL_SLACKS_POS', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, folder_path, "EL_Slacks Demand")

    # Save heat demand
    save_detailed(MODEL_INPUT.HEAT_DEMAND.TIMESERIES, names(INPUT.IN_HEAT_DEMAND_TIMESERIES)[2:end],
        TMP_CALC_TIME, folder_path, "TH_Heat Demand")
    # Save heat generation of CHP units
    save_detailed(RESULTS.GEN_HEAT_CHP', INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.unit_name],
        TMP_CALC_TIME, folder_path, "TH_Gen CHP")
    # Save generation of heat only units
    save_detailed(RESULTS.GEN_HEAT_ONLY', INPUT.IN_HEAT_ONLY[:,COLUMN_DEFINITIONS.COLUMNS_HEAT_ONLY.unit_name],
        TMP_CALC_TIME, folder_path, "TH_Gen Heat Only")
    # Save heat region slacks
    save_detailed(RESULTS.SLACK_Q', INPUT.IN_HEAT_NODE[:,COLUMN_DEFINITIONS.COLUMNS_HEATREGION.heatregion_id],
        TMP_CALC_TIME, folder_path, "TH_Heat Slacks")

    # Save hydrogen generation of PTG units
    save_detailed((RESULTS.GEN_PTG .* MODEL_INPUT.PTG.ETA_P_MAX)', INPUT.IN_PTG[:,COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.unit_name],
        TMP_CALC_TIME, folder_path, "H2_Gen PTG")

    # Save prices
    save_detailed(RESULTS.PRICES', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, folder_path, "EL_LMPs")

    #-- Save AC branches' loading
    # Get columns with prefixes
    el_branch_columns = vcat("Transformer_".*INPUT.IN_TRAFOS[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name],
        "Line_".*INPUT.IN_LINES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name])
    # Save line loading
    save_detailed(RESULTS.FLOW_EL_BRANCH', el_branch_columns,
        TMP_CALC_TIME, folder_path, "EL_AC Branches (Lines & Transformers)")

    # Get set of AC lines
    SET_AC_LINES = MODEL_INPUT.SETS.BRANCHES[size(INPUT.IN_TRAFOS,1)+1:end]
    # Save line loading
    save_detailed(RESULTS.FLOW_EL_BRANCH[SET_AC_LINES,:]', INPUT.IN_LINES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name],
        TMP_CALC_TIME, folder_path, "EL_AC Lines")

    # Get set of AC lines
    SET_AC_TRAFOS = MODEL_INPUT.SETS.BRANCHES[1:size(INPUT.IN_TRAFOS,1)]
    # Save line loading
    save_detailed(RESULTS.FLOW_EL_BRANCH[SET_AC_TRAFOS,:]', INPUT.IN_TRAFOS[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name],
        TMP_CALC_TIME, folder_path, "EL_AC Transformers")

    # Get set of AC lines
    SET_REL_AC_LINES = MODEL_INPUT.SETS.BRANCHES[size(INPUT.IN_TRAFOS,1)+1:end]
    # Get absolute values of relative power flows, i.e. normalized using the rated power
    ABS_FLOW_EL_LINE = broadcast(abs, RESULTS.FLOW_EL_BRANCH[SET_REL_AC_LINES,:] ./ MODEL_INPUT.BRANCHES.PMAX[SET_REL_AC_LINES])

    # Save AC line loading
    save_detailed(ABS_FLOW_EL_LINE', el_branch_columns[SET_REL_AC_LINES],
        TMP_CALC_TIME, folder_path, "EL_AC Lines (rel.)")

    #-- Save E-Mob
    # Save E-Mob residential
    save_detailed(MODEL_INPUT.EMOB_RESIDENTIAL.TIMESERIES, names(INPUT.IN_EMOB_RESIDENTIAL_TIMESERIES)[2:end],
        TMP_CALC_TIME, folder_path, "EL_EMOB_residential")
    # Save E-Mob
    save_detailed(MODEL_INPUT.EMOB_COMMERCIAL.TIMESERIES, names(INPUT.IN_EMOB_COMMERCIAL_TIMESERIES)[2:end],
        TMP_CALC_TIME, folder_path, "EL_EMOB_comercial")
    # Save E-Mob
    save_detailed(MODEL_INPUT.EMOB_QUICK.TIMESERIES, names(INPUT.IN_EMOB_QUICK_TIMESERIES)[2:end],
        TMP_CALC_TIME, folder_path, "EL_EMOB_quickcharge")

end


function preprocess_output(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    # Get quantiles of brnach loading
    QUANT_ABS_FLOW_EL_LINE = quantiles_line_loading(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS)
    # Define reults matrices and names
    results = Dict(
        1 => MODEL_INPUT.LOADS.TIMESERIES',
        2 => RESULTS.GEN_HEAT_PUMP,
        3 => RESULTS.GEN_HEAT_PTH,
        4 => RESULTS.GEN_PTG,
        5 => MODEL_INPUT.RENEWABLES.TIMESERIES[:,1:size(INPUT.IN_PV,1)]',
        6 => MODEL_INPUT.RENEWABLES.TIMESERIES[:,size(INPUT.IN_PV,1)+1:end]',
        7 => MODEL_INPUT.HY_ROR.TIMESERIES',
        8 => RESULTS.RENEWABLE_FEED_IN,
        9 => RESULTS.GEN_EL_BIO,
        10 => RESULTS.GEN_EL_ONLY,
        11 => RESULTS.GEN_EL_CHP,
        12 => RESULTS.STORAGES_CHARGING,
        13 => RESULTS.STORAGES_DISCHARGING,
        14 =>  RESULTS.HVDC_FLOW,
        15 => RESULTS.GEN_EL_SLACKS_NEG,
        16 => RESULTS.GEN_EL_SLACKS_POS,
        17 => MODEL_INPUT.HEAT_DEMAND.TIMESERIES',
        18 => MODEL_INPUT.CONST_HEAT_LOAD.Q_MAX,
        19 => RESULTS.GEN_HEAT_CHP,
        20 => RESULTS.GEN_HEAT_ONLY,
        21 => RESULTS.SLACK_Q,
        22 => RESULTS.GEN_PTG' * MODEL_INPUT.PTG.ETA_P_MAX,
        23 => RESULTS.PRICES,
        24 => RESULTS.FLOW_EL_BRANCH,
        25 => QUANT_ABS_FLOW_EL_LINE'
        )

    results_names = Dict(
        1 => "El. Load (conv.)",
        2 => "El. Load (HP)",
        3 => "El. Load (PTH)",
        4 => "El. Load (PTG)",
        5 => "PV",
        6 => "Wind",
        7 => "Hydro ROR",
        8 => "RE (curtailed)",
        9 => "Bio",
        10 => "El. Gen (Gas)",
        11 => "El. Gen (CHP)",
        12 => "El. Storages (Charge)",
        13 => "El. Storages (Discharge)",
        14 => "HVDC Lines",
        15 => "El. Slack (neg.)",
        16 => "El. Slack (pos.)",
        17 => "Th. Demand",
        18 => "Th. Const. Demand",
        19 => "Th. Gen. (CHP)",
        20 => "Th. Gen. (Heat Only)",
        21 => "Th. Slack",
        22 => "H2 Gen (PTG)",
        23 => "Prices",
        24 => "AC Lines",
        25 => "Quantiles AC Lines"
        )
    return results, results_names
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

function fill_sheet(sheet, arr)
    for ind in CartesianIndices(arr)
        XLSX.setdata!(sheet, XLSX.CellRef(ind[1], ind[2]), arr[ind])
    end
end

function save_variable2excel(file_name, var, sheet_name)
    XLSX.openxlsx(file_name*".xlsx", mode="w") do xf
        s  = XLSX.addsheet!(xf,sheet_name)
        fill_sheet(s,var)
    end
end

function save_dataframe2excel(file_name, df, sheet_name)
    XLSX.openxlsx(file_name*".xlsx", mode="w") do xf
        XLSX.addsheet!(xf,sheet_name)
        sheet = xf[sheet_name]
        for r in 1:size(df,1), c in 1:size(df,2)
          sheet[XLSX.CellRef(r , c )] = df[r,c]
        end
    end
end

function save_dataframe2excel_add(file_name, df, sheet_name)
    XLSX.openxlsx(file_name*".xlsx", mode="rw") do xf
        XLSX.addsheet!(xf,sheet_name)
        sheet = xf[sheet_name]
        for r in 1:size(df,1), c in 1:size(df,2)
          sheet[XLSX.CellRef(r , c )] = df[r,c]
        end
    end
end



#--- get version function
function get_version_info(in_names::Vector{String},in_placements::Vector{UnitRange{Int64}}, in_format_version::Regex)
    # this function returns a version info (if there is one) or otherwise "no_version"

    #copy input names
    tmp_in_names = copy(in_names)
    #get index of end of format of technology/unit
    placements_end = last.(in_placements)
    #get lenght of complete names (to  know the last index in the string)
    names_length = length.(tmp_in_names)
    #cut off the tail-element of the string after the end of format of technology/unit
    cut_names_old = getindex.(tmp_in_names, range.(placements_end,names_length, step=1))
    #find actual version in the cut off tail (string)
    placements_version_info = findfirst.(in_format_version,cut_names_old)
    #check if there is no version ending end return "no_version" or return actual version endings if there are
    if (size(unique(placements_version_info),1) == 1) & (unique(placements_version_info)[1] == nothing)
        unique_version_info = ["no_version",]
    else
        begin_placements_version_info = first.(placements_version_info)
        version_info = getindex.(cut_names_old,range.(begin_placements_version_info,length.(cut_names_old),step=1))
        unique_version_info = unique(version_info)
    end
    #place a warning if there are multiple versions on the same technology
    if size(unique(unique_version_info),1) > 1
        println("WARNING: There is more than one version for the same Technology")
    end

    #return version info
    return unique_version_info
end

#--- rename function
function rename_units_for_output(in_DF::DataFrame,in_format::Regex,in_format_version::Regex)
    # this function renames the units in the column names of a DF "in_DF" in the format of "in_format"
    #if there ist a version ending in the given format the version info is returned as second value

    #copy input DF to not modify the original DF
    tmp_DF = copy(in_DF)
    #get names from DF
    #check if input DF has Timesteps and print exception if not
    if names(tmp_DF)[1] == "Timestep"
        in_names = names(tmp_DF)[2:end]
    else
        in_names = names(tmp_DF)
        println("WARNING: Input DF has no column Timesteps")
    end

    #throw error if column names are not unique
    if size(in_names,1) !== size(unique(in_names),1)
        error("ERROR: Column names of DF are not unique")
    end
    #find the placement of the given format "in_format" in the old names
    placements = findfirst.(in_format,in_names)
    #throw error if given in_format format is not present in strings
    if nothing in placements
        error("ERROR: Could not find given format in at least one of the DF column names")
    end
    #get version info
    version_info = get_version_info(in_names, placements, in_format_version)

    #workaround for WIND_OFF (Since WIND_OFF is far away from WO and close to WIND (string wise) and also handled at the same time as WIND)
    tmp_in_names = replace.(in_names, "WIND_OFF" => "WO")
    out_names = getindex.(tmp_in_names,placements)
    #if "WIND_OFF" is named as "WO" the previous 2 lines can be reformulated as: out_names = getindex.(in_names,placements)

    #throw error if column names are not unique
    if size(in_names,1)!==size(unique(out_names),1)
        error("ERROR: Cant make column names of DF for output unique - this might happen e.g. if multiple technologies are put together for an output and get the same name")
    end
    #create dictionary for renaming units (column names) to the new format
    rename_dict = Dict(in_names .=> out_names)
    #rename units (column names) and save result in out_DF
    out_DF = rename(tmp_DF,rename_dict)
    #return new DF
    return out_DF,version_info
end


#--- function to create DFs for online model output
function create_DF_for_online_model(IN_DF, col_names, time)
    # Get column names
    tmp_cols = broadcast(string,vcat("Timestep",col_names))
    # Get data as matrix
    tmp_mat = hcat(time, Matrix{Float64}(IN_DF))
    # Create DataFrame
    tmp_df = DataFrame(tmp_mat, tmp_cols)
    return tmp_df
end

#--- function to save DF as file_name.csv in given (relative) folder_path
function save_DF_as_csv_in_path(IN_DF,file_name::String,folder_path::String)
    # Save data to csv
    CSV.write(folder_path*"\\"*file_name*".csv", IN_DF)
end

#--- function to save version info to exising .txt or create new if not existing
function save_version_to_existing_txt_in_path(version_info::Vector{String}, file_name::String, folder_path::String)
    print_version_info = join(version_info, " ")
    open(folder_path * "\\" * "Version_Info.txt", "a") do io
       write(io, file_name * " " .* print_version_info * "\n")
    end
end

#--- function to save data for online model
function save_data_for_online_model(IN_DF, col_names, time, file_name, format, version_format, folder_path)
    tmp_DF = create_DF_for_online_model(IN_DF, col_names, time)
    tmp_DF, tmp_version = rename_units_for_output(tmp_DF,format,version_format)
    save_DF_as_csv_in_path(tmp_DF,file_name,folder_path)
    save_version_to_existing_txt_in_path(tmp_version,file_name,folder_path)
end

#--- output for online model
function save_detailed_results_for_online_model(RESULTS, MODEL_INPUT, INPUT, PARAMETER_SETTINGS, COLUMN_DEFINITIONS)
    #get index of timeseries input as temporary calc time
    TMP_CALC_TIME = collect(1:length(PARAMETER_SETTINGS.CALC_TIME))
    #get formating informations
    format = parameters_names_and_settings.format_definitions()

    # Define folder name
    folder_path = abspath(joinpath(pwd(), "../results/")) * Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * " - ResultsForOnlineModel"
    # Create folder if not exists
    mkdir(folder_path)

    # Save electrical load
    save_data_for_online_model(MODEL_INPUT.LOADS.TIMESERIES, names(INPUT.IN_LOAD_TIMESERIES)[2:end],
        TMP_CALC_TIME, "EL_Load(conv.)",
        format.loads,format.version,folder_path)
    # Save commercial BEV's electrical load
    save_data_for_online_model(MODEL_INPUT.BEV.COMMERCIAL_TIMESERIES, names(INPUT.IN_BEV_COMMERCIAL_ELEC_TIMESERIES)[2:end],
        TMP_CALC_TIME, "Load (BEV commercial)", format.loads,format.version,folder_path)
    # Save quickcharging BEV's electrical load
    save_data_for_online_model(MODEL_INPUT.BEV.QUICKCHARGE_TIMESERIES, names(INPUT.IN_BEV_QUICKCHARGE_ELEC_TIMESERIES)[2:end],
        TMP_CALC_TIME, "Load (BEV quickcharge)", format.loads,format.version,folder_path)
    # Save residential BEV's electrical load
    save_data_for_online_model(MODEL_INPUT.BEV.RESIDENTIAL_TIMESERIES, names(INPUT.IN_BEV_RESIDENTIAL_ELEC_TIMESERIES)[2:end],
        TMP_CALC_TIME, "Load (BEV residential)", format.loads,format.version,folder_path)
    # Save heat pump's load
    save_data_for_online_model(RESULTS.GEN_HEAT_PUMP', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_Load(HP)",
        format.loads,format.version,folder_path)
    # Save heat PTH's load
    save_data_for_online_model(RESULTS.GEN_HEAT_PTH', INPUT.IN_PTH_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_PTH.unit_name],
        TMP_CALC_TIME, "EL_Load(PTH)",
        format.loads,format.version,folder_path)
    # Save heat PTG's load
    save_data_for_online_model(RESULTS.GEN_PTG', INPUT.IN_PTG[:,COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.unit_name],
        TMP_CALC_TIME, "EL_Load(PTG)",
        format.loads,format.version,folder_path)
    # Save PV generation
    save_data_for_online_model(MODEL_INPUT.RENEWABLES.TIMESERIES[:,1:size(INPUT.IN_PV,1)], names(INPUT.IN_RENEWABLE_TIMESERIES)[2:size(INPUT.IN_PV,1)+1],
        TMP_CALC_TIME, "EL_PV",
        format.renewables,format.version,folder_path)
    # Save Wind generation
    save_data_for_online_model(MODEL_INPUT.RENEWABLES.TIMESERIES[:,size(INPUT.IN_PV,1)+1:end], names(INPUT.IN_RENEWABLE_TIMESERIES)[2+size(INPUT.IN_PV,1):end],
        TMP_CALC_TIME, "EL_Wind",
        format.renewables,format.version,folder_path)
    # Save hydro ron-of-the-river generation
    save_data_for_online_model(MODEL_INPUT.HY_ROR.TIMESERIES, names(INPUT.IN_HY_ROR_TIMESERIES)[2:end],
        TMP_CALC_TIME, "EL_Hydro-ROR",
        format.renewables,format.version,folder_path)
    # Save curtailed renewable feed-in
    save_data_for_online_model(RESULTS.RENEWABLE_FEED_IN', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_Renewable_curtailed(PV,Wind,Hydro-ROR)",
        format.nodes,format.version,folder_path)
    # Save biomass/biogas generation
    save_data_for_online_model(RESULTS.GEN_EL_BIO', INPUT.IN_BIO[:,COLUMN_DEFINITIONS.COLUMNS_BIO.unit_name],
        TMP_CALC_TIME, "EL_Biomass,Biogas",
        format.renewables,format.version,folder_path)
    # Save electricty only generation
    save_data_for_online_model(RESULTS.GEN_EL_ONLY', INPUT.IN_CONVENTIONAL[:,COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.unit_name],
        TMP_CALC_TIME, "EL_Gen_Elec_Only",
        format.conventionals,format.version,folder_path)
    # Save CHP generation
    save_data_for_online_model(RESULTS.GEN_EL_CHP', INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.unit_name],
        TMP_CALC_TIME, "EL_Gen_CHP",
        format.conventionals,format.version,folder_path)
    # Save pumped storages charging
    save_data_for_online_model(RESULTS.STORAGES_CHARGING', INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.unit_name],
        TMP_CALC_TIME, "EL_Pumped-storage_charging",
        format.renewables,format.version,folder_path)
    # Save pumped storages discharging
    save_data_for_online_model(RESULTS.STORAGES_DISCHARGING', INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.unit_name],
        TMP_CALC_TIME, "EL_Pumped-storage_discharging",
        format.renewables,format.version,folder_path)
    # Save emob residential
    save_data_for_online_model(RESULTS.LOAD_EMOB_RESIDENTIAL', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_Emob_residential",
        format.loads, format.version,folder_path)
    save_data_for_online_model(RESULTS.LOAD_EMOB_COMMERCIAL', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_Emob_commercial",
        format.loads, format.version,folder_path)
    save_data_for_online_model(RESULTS.LOAD_EMOB_QUICKCHARGE', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_Emob_quickcharge",
        format.loads, format.version,folder_path)


    # ==== SLACKS =====
    # Save HVDC flow
    save_detailed(RESULTS.HVDC_FLOW', INPUT.IN_HVDC[:,COLUMN_DEFINITIONS.COLUMNS_HVDC.unit_name],
        TMP_CALC_TIME, folder_path, "EL_HVDC Flow")
    # Save slack generation
    save_data_for_online_model(RESULTS.GEN_EL_SLACKS_NEG', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_Slacks_Generation",
        format.nodes,format.version,folder_path)
    # Save slack demand
    save_data_for_online_model(RESULTS.GEN_EL_SLACKS_POS', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_Slacks_Demand",
        format.nodes,format.version,folder_path)

    # Save prices
    save_data_for_online_model(RESULTS.PRICES', INPUT.IN_NODE[:,COLUMN_DEFINITIONS.COLUMNS_NODE.node_name],
        TMP_CALC_TIME, "EL_LMPs",
        format.nodes,format.version,folder_path)

    # Get set of AC lines
    SET_AC_LINES = MODEL_INPUT.SETS.BRANCHES[size(INPUT.IN_TRAFOS,1)+1:end]
    # Save line loading
    save_data_for_online_model(RESULTS.FLOW_EL_BRANCH[SET_AC_LINES,:]', INPUT.IN_LINES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name],
        TMP_CALC_TIME, "EL_AC_Lines",
        format.lines,format.version,folder_path)

    # Get set of Trafos
    SET_AC_TRAFOS = MODEL_INPUT.SETS.BRANCHES[1:size(INPUT.IN_TRAFOS,1)]
    # Save line loading
    save_data_for_online_model(RESULTS.FLOW_EL_BRANCH[SET_AC_TRAFOS,:]', INPUT.IN_TRAFOS[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name],
        TMP_CALC_TIME, "EL_AC_Transformers",
        format.trafos,format.version,folder_path)

    # Get set of AC lines
    SET_REL_AC_LINES = MODEL_INPUT.SETS.BRANCHES[size(INPUT.IN_TRAFOS,1)+1:end]
    # Get absolute values of relative power flows, i.e. normalized using the rated power
    ABS_FLOW_EL_LINE = broadcast(abs, RESULTS.FLOW_EL_BRANCH[SET_REL_AC_LINES,:] ./ MODEL_INPUT.BRANCHES.PMAX[SET_REL_AC_LINES])

    # Save AC line loading
    el_branch_columns = vcat("Transformer_".*INPUT.IN_TRAFOS[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name],
        "Line_".*INPUT.IN_LINES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.name])
    
    save_data_for_online_model(ABS_FLOW_EL_LINE', el_branch_columns[SET_REL_AC_LINES],
        TMP_CALC_TIME, "EL_AC_Lines(rel.)",
        format.lines,format.version,folder_path)

    return nothing

end
#---
function save_results_for_AoTRF(MODEL_INPUT,RESULTS)
    RELEVANT_INPUT=(
        P_MAX_HP=MODEL_INPUT.HEATPUMPS.P_MAX_HP,
        COP_TIMESERIES=MODEL_INPUT.HEATPUMPS.COP_TIMESERIES,
        P_MAX_HE=MODEL_INPUT.HEATPUMPS.P_MAX_HE,
        Q_MAX=MODEL_INPUT.HEATPUMPS.Q_MAX,
        NU_SELF_DIS=MODEL_INPUT.HEATPUMPS.NU_SELF_DIS,
        ELEC_TIMESERIES=MODEL_INPUT.HEATPUMPS.ELEC_TIMESERIES,
        HEAT_TIMESERIES=MODEL_INPUT.HEATPUMPS.HEAT_TIMESERIES,
    )

    #save as xlsx
    save_results_as_excel_ffs_2(RESULTS, "Results_AoTRF")
    save_results_as_excel_ffs_2(RELEVANT_INPUT, "RELEVANT_INPUT_AoTRF")

end

#---
function save_as_jld2(DATA, NAME)
    str_name = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS") * NAME * ".jld2"
    save_object(str_name,DATA)
end


function save_multiple_jld2(DATA, NAME)
    timestring = Dates.format(now(), "yyyy_mm_dd_HH_MM_SS")
    for ct_file in NAME
        idx = findall(x -> x == ct_file, NAME)
        file_path = abspath(joinpath(pwd(), "../results/") * timestring * " - " * ct_file * ".jld2")
        save_object(file_path,DATA[idx])
    end
    return nothing
end


end

