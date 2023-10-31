module input_functions

using XLSX, CSV, DataFrames, Setfield, SparseArrays, InteractiveUtils 

include("basic_functions.jl")
include("load_data.jl")

#--- Preprocess false datatypes of strings
function change_datatype_strings(df)
    try
        columns = names(df)
        strings = subtypes(InlineString)
        for col in columns
            if eltype(df[!,col]) in strings
                df[!,col] = convert.(String,df[!,col])
            end
        end
    finally
        return copy(df)
    end
end


#--- Preprocess renewable assignment data
function preprocess_renewable_data(data_assignment)
    # Initialize new DataFrame holding assignment data
    IN_RENEWABLE = copy(change_datatype_strings(data_assignment.elec_pv))
    # Append rows below each other
    DataFrames.append!(IN_RENEWABLE, change_datatype_strings(data_assignment.elec_wind))

    IN_PV = data_assignment.elec_pv
    IN_WIND = data_assignment.elec_wind
    return IN_RENEWABLE, IN_PV, IN_WIND
end

#--- Preprocess renewable timeseries data
function preprocess_renewable_timeseries_data(data_timeseries, CALC_TIME)
    # Initialize matrix holding renewable timeseries data
    IN_RENEWABLE_TIMESERIES = copy(data_timeseries.elec_pv)
    # Get column names of wind timeseries
    tmp_names = names(data_timeseries.elec_wind)
    # Append columns of wind timeseries data
    [IN_RENEWABLE_TIMESERIES[!,col] = data_timeseries.elec_wind[:,col] for col in tmp_names[2:end]]
    # Limit renewable timeseries to considered time of calculation
    IN_RENEWABLE_TIMESERIES = IN_RENEWABLE_TIMESERIES[CALC_TIME,:]
    return IN_RENEWABLE_TIMESERIES
end

#--- Preprocess nodal assignment and unit parameters
function join_parameter_assignment(DF_PARAMETER,DF_ASSIGNMENT)
    # assigns the corrseponding node from DF_ASSIGNMENT to the unit in DF_PARAMETER and returns a new Dataframe with parameers and nodes
    # creating copy of Input DataFrame DF_PARAMETER & sorting it
    tmp_sorted_df = copy(DF_PARAMETER)
    sort!(tmp_sorted_df, :unit_name)
    # Error handling if Input DataFrame DF_PARAMETER is not sorted
    if tmp_sorted_df.unit_name != DF_PARAMETER.unit_name
        error("Input DataFrame is not sorted! The Input DataFrame needs to be sorted by unit_name to sort it again after the join.")
    end
    DF_JOINED = leftjoin(DF_PARAMETER,DataFrame(unit_id=DF_ASSIGNMENT.unit_id, node_id=DF_ASSIGNMENT.node_id), on = :unit_id)
    # The Input DataFrame needs to be sorted by unit_name to sort it again after the join
    sort!(DF_JOINED, :unit_name)
    return DF_JOINED
end
#--- Preprocess nodal assignment and unit parameters
function join_parameter_assignment_heat(DF_PARAMETER,DF_ASSIGNMENT)
    # assigns the corrseponding node from DF_ASSIGNMENT to the unit in DF_PARAMETER and returns a new Dataframe with parameers and nodes
    # creating copy of Input DataFrame DF_PARAMETER & sorting it
    tmp_sorted_df = copy(DF_PARAMETER)
    sort!(tmp_sorted_df, :unit_name)
    # Error handling if Input DataFrame DF_PARAMETER is not sorted
    if tmp_sorted_df.unit_name != DF_PARAMETER.unit_name
        error("Input DataFrame is not sorted! The Input DataFrame needs to be sorted by unit_name to sort it again after the join.")
    end
    DF_JOINED = leftjoin(DF_PARAMETER, DataFrame(unit_id=DF_ASSIGNMENT.unit_id, heat_node_id=DF_ASSIGNMENT.heat_node_id), on = :unit_id)
    # The Input DataFrame needs to be sorted by unit_name to sort it again after the join
    sort!(DF_JOINED, :unit_name)
    return DF_JOINED
end

function preprocess_heatpump_timeseries_data(data_timeseries, CALC_TIME)
    # Initialize matrix holding heat pump electricity demand timeseries data
    IN_HEATPUMP_TIMESERIES = copy(data_timeseries.elec_hp_heat_demand)
    # Get COP timeseries of heat pumps
    TMP_HEATPUMP_COP_TIMESERIES = copy(data_timeseries.elec_hp_cop)
    # Get column names of COP timeseries
    tmp_names_cop = names(data_timeseries.elec_hp_cop)
    # Get column names of heat demand timeseries
    tmp_names_heat_demand = names(data_timeseries.elec_hp_heat_demand)
    #print(tmp_names_heat_demand)
    # Replace COP in name of elec heat pump
    # tmp_names_cop = [replace(col, "COP_" => "DEMAND_") for col in tmp_names_cop]
    tmp_names_cop = [replace(col, "COP_" => "THERMAL_") for col in tmp_names_cop]
    # Replace DEMAND in name of elec heat pump
    tmp_names_heat_demand = [replace(col, "DEMAND_" => "") for col in tmp_names_heat_demand]
    # Rename columns of heat demand timeseries
    DataFrames.rename!(IN_HEATPUMP_TIMESERIES, Symbol.(tmp_names_heat_demand))
    # Rename columns of COP timeseries
    DataFrames.rename!(TMP_HEATPUMP_COP_TIMESERIES, Symbol.(tmp_names_cop))
    # Get heat pump electricty demand timeseries
    for col in tmp_names_heat_demand[2:end]
        IN_HEATPUMP_TIMESERIES[:, col] = IN_HEATPUMP_TIMESERIES[:, col] ./ TMP_HEATPUMP_COP_TIMESERIES[:, col]
    end
    return IN_HEATPUMP_TIMESERIES
end


function rename_and_sort_timeseries_DFa_names_by_DFb_column(DFa, DFb_column)
    tmp_DFa = copy(DFa)
    tmp_DFb_column=copy(DFb_column)
    new_names=replace.(names(tmp_DFa),"_DEMAND"=>"")
    new_names=replace.(new_names,"_THERMAL"=>"")
    DataFrames.rename!(tmp_DFa,Symbol.(new_names))
    DFreturn = tmp_DFa[!,Symbol.(tmp_DFb_column)]
    insertcols!(DFreturn,1,:timestep=>tmp_DFa.timestep)
    return DFreturn
end

function load_input_data(COLUMN_DEFINITIONS, PARAMETER_SETTINGS)
    # Load data from csv files
    data_timeseries, data_unit_parameters, data_assignment, data_elec_network, data_heat_network = load_data.load_data_from_tmp_folder(PARAMETER_SETTINGS)
    # Prepare branch parameters for load flow
    IN_BRANCHES, IN_TRAFOS, IN_LINES  = basic_functions.preprocess_ac_lines_and_transformers(data_elec_network,
        COLUMN_DEFINITIONS.COLUMNS_TRANSFORMER, COLUMN_DEFINITIONS.COLUMNS_LINE,
        COLUMN_DEFINITIONS.COLUMNS_NODE, PARAMETER_SETTINGS.REF_VOLTAGE)
    # overwrite low voltage lines with arbitrarily high line Limit
    #idx_low_voltage_lovel = (COLUMN_DEFINITIONS.COLUMNS_LINE.voltage_level < PARAMETER_SETTINGS.LINE_RESTRICTIONS)
    # get hvdc table
    IN_HVDC = data_elec_network.hvdc_lines
    # Get conventional's nodal assignment and unit parameters in one dataframe
    IN_CONVENTIONAL = join_parameter_assignment(
        basic_functions.get_technical_parameters(data_unit_parameters.elec_only_units,
        COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL),data_assignment.elec_only_units)
    # Get storage's nodal assignment and unit parameters in one dataframe
    IN_STORAGE = join_parameter_assignment(data_unit_parameters.el_storage_units, data_assignment.el_storage_units)
    # Filter only for considered storage units in relevant voltage levels
    IN_STORAGE = leftjoin(IN_STORAGE, data_elec_network.nodes[:,[:node_id, :voltage_id]], on=:node_id )
    idx_relev_voltagelevel = map(in(Set(PARAMETER_SETTINGS.FLEXIBILITY_STORAGE_VOLTAGE_LEVELS)), IN_STORAGE.voltage_id)
    idx_relev_battery = occursin.("_PS_", IN_STORAGE.unit_name)
    idx_relev = (idx_relev_voltagelevel .| idx_relev_battery)
    println("Only consider STOR in voltage levels " * join(PARAMETER_SETTINGS.FLEXIBILITY_STORAGE_VOLTAGE_LEVELS, ", ") * " : From all " * string(size(idx_relev,1)) * " BS-STOR units will remain " * string(sum(idx_relev)) *  " STOR units.")
    IN_STORAGE = IN_STORAGE[idx_relev,:]
    # Set self discharging ratio to 0.1 %
    IN_STORAGE.self_dis_ratio = ones(size(IN_STORAGE,1))*0.001
    
    # Get meta data of PV and Wind in one dataframe
    IN_RENEWABLE, IN_PV, IN_WIND = preprocess_renewable_data(data_assignment)
    # Get PV and Wind timeseries data in one dataframe
    IN_RENEWABLE_TIMESERIES = preprocess_renewable_timeseries_data(data_timeseries,
        PARAMETER_SETTINGS.CALC_TIME)
    # Get meta data of heat pumps in one dataframe
    IN_HEATPUMP = join_parameter_assignment(data_unit_parameters.elec_hp,
        data_assignment.elec_hp)
    # Emob
    IN_HEATPUMP = join_parameter_assignment(data_unit_parameters.elec_hp, data_assignment.elec_hp)


    # Get heat pump's electricty demand timeseries in one dataframe
    IN_HEATPUMP_ELEC_TIMESERIES = preprocess_heatpump_timeseries_data(data_timeseries, PARAMETER_SETTINGS.CALC_TIME)
    IN_HEATPUMP_ELEC_TIMESERIES = rename_and_sort_timeseries_DFa_names_by_DFb_column(IN_HEATPUMP_ELEC_TIMESERIES, IN_HEATPUMP.unit_name)
    # Get heat pump's heat demand timeseries in one dataframe
    IN_HEATPUMP_HEAT_TIMESERIES = copy(data_timeseries.elec_hp_heat_demand)
    IN_HEATPUMP_HEAT_TIMESERIES = rename_and_sort_timeseries_DFa_names_by_DFb_column(IN_HEATPUMP_HEAT_TIMESERIES, IN_HEATPUMP.unit_name)
    # Get heat pump's COP timeseries in one dataframe
    IN_HEATPUMP_COP_TIMESERIES = copy(data_timeseries.elec_hp_cop)
    # Get CHP unit's nodal assignment and unit parameters in one dataframe
    IN_COMBINEDHEATPOWER_ELEC = join_parameter_assignment(basic_functions.get_technical_parameters(
        data_unit_parameters.chp_units, COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER),data_assignment.chp_units_elec)
    sort!(IN_COMBINEDHEATPOWER_ELEC, COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.unit_name)
    # Get CHP unit's heat region assignment as dataframe
    IN_COMBINEDHEATPOWER_HEAT = data_assignment.chp_units_heat
    sort!(IN_COMBINEDHEATPOWER_HEAT, COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.unit_name)
    # Power-to-gas
    IN_PTG = join_parameter_assignment(basic_functions.get_technical_parameters(
        data_unit_parameters.ptg_hydrogen_units, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS),data_assignment.ptg_hydrogen_units)

    INPUT = (
        IN_NODE = data_elec_network.nodes,
        IN_BRANCHES,
        IN_LINES,
        IN_TRAFOS,
        IN_HVDC,
        IN_PTG,
        IN_CONVENTIONAL,
        IN_RENEWABLE,
        IN_PV,
        IN_WIND,
        IN_RENEWABLE_TIMESERIES,
        IN_LOAD = data_assignment.elec_demand,
        IN_LOAD_TIMESERIES = data_timeseries.elec_demand[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_STORAGE,
        IN_COMBINEDHEATPOWER_ELEC,
        IN_COMBINEDHEATPOWER_HEAT,
        IN_BIO = data_assignment.elec_bio,
        IN_BIO_TIMESERIES = data_timeseries.elec_bio[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_HY_ROR = data_assignment.elec_hy_ror,
        IN_HY_ROR_TIMESERIES = data_timeseries.elec_hy_ror[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_EMOB_RESIDENTIAL = data_assignment.elec_bev_residential,
        IN_EMOB_COMMERCIAL = data_assignment.elec_bev_commercial,
        IN_EMOB_QUICK = data_assignment.elec_bev_quick,
        IN_EMOB_RESIDENTIAL_TIMESERIES = data_timeseries.elec_bev_residential[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_EMOB_RESIDENTIAL_TIMESERIES_PMAX = data_timeseries.elec_bev_residential_pmax[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_EMOB_COMMERCIAL_TIMESERIES = data_timeseries.elec_bev_commercial[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_EMOB_QUICK_TIMESERIES = data_timeseries.elec_bev_quick[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_HEATPUMP,
        IN_HEATPUMP_ELEC_TIMESERIES = IN_HEATPUMP_ELEC_TIMESERIES[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_HEATPUMP_HEAT_TIMESERIES = IN_HEATPUMP_HEAT_TIMESERIES[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_HEATPUMP_COP_TIMESERIES = IN_HEATPUMP_COP_TIMESERIES[PARAMETER_SETTINGS.CALC_TIME,:],
        IN_HEAT_NODE = data_heat_network.nodes,
        IN_PTH_ELEC = join_parameter_assignment(data_unit_parameters.heat_pth, data_assignment.elec_pth),
        IN_PTH_HEAT = join_parameter_assignment_heat(data_unit_parameters.heat_pth, data_assignment.heat_pth),
        IN_HEAT_ONLY = join_parameter_assignment_heat(data_unit_parameters.heat_only, data_assignment.heat_only),
        IN_CONST_HEAT_LOAD = join_parameter_assignment_heat(data_unit_parameters.const_heat_load, data_assignment.const_heat_load),
        IN_HEAT_DEMAND = data_assignment.heat_demand,
        IN_HEAT_DEMAND_TIMESERIES = data_timeseries.heat_demand[PARAMETER_SETTINGS.CALC_TIME,:],
    )
    return INPUT
end

function fill_technology_structs(INPUT, COLUMN_DEFINITIONS, PARAMETER_SETTINGS)

    println("----------")
    #--- Create assignment matrices for mapping conventionals, renewables and laod to nodes
    #-- Power-to-gas assignment matrix (Dimension: nodes x powerplants)
    A_PTG = basic_functions.create_assignment_matrix(INPUT.IN_NODE, INPUT.IN_PTG, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.node_id, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- Conventional A_G assignment matrix (Dimension: nodes x powerplants)
    A_G = basic_functions.create_assignment_matrix(INPUT.IN_NODE, INPUT.IN_CONVENTIONAL, COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.node_id, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- Renewable A_REN assignment matrix (Dimension: nodes x powerplants)
    A_REN = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_RENEWABLE, INPUT.IN_RENEWABLE_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_RENEWABLE.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- Load A_L assignment matrix (Dimension: nodes x loads)
    A_L = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_LOAD, INPUT.IN_LOAD_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_LOAD.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- storage A_L assignment matrix (Dimension: nodes x storages)
    A_STOR = basic_functions.create_assignment_matrix(INPUT.IN_NODE, INPUT.IN_STORAGE, COLUMN_DEFINITIONS.COLUMNS_STORAGE.node_id, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- Combined heat power plants (CHP) A_CHP assignment matrix (Dimension: nodes x CHP)
    A_CHP_EL = basic_functions.create_assignment_matrix(INPUT.IN_NODE, INPUT.IN_COMBINEDHEATPOWER_ELEC, COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.node_id, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- HEATload A_CHP_H assignment matrix (Dimension: nodes x heatloads)
    A_CHP_H = basic_functions.create_assignment_matrix(INPUT.IN_HEAT_NODE, INPUT.IN_COMBINEDHEATPOWER_HEAT, COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.heatregion_id, COLUMN_DEFINITIONS.COLUMNS_HEATREGION.heatregion_id)
    #-- HVDC assignment matrix (Dimension: HVDC lines x nodes)
    A_HVDC, A_HVDC_LOSS = basic_functions.create_assignment_matrix_HVDC(INPUT.IN_HVDC, INPUT.IN_NODE, COLUMN_DEFINITIONS.COLUMNS_NODE, COLUMN_DEFINITIONS.COLUMNS_HVDC)
    #-- Biomass/gas assignment matrix (Dimension: nodes x powerplants)
    A_BIO = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_BIO, INPUT.IN_BIO_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_BIO.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    # basic_functions.create_assignment_matrix(INPUT.IN_NODE, INPUT.IN_BIO, COLUMN_DEFINITIONS.COLUMNS_BIO.node_id, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- Hdyro run-of-river assignment matrix (Dimension: nodes x powerplants)
    A_HY_ROR = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_HY_ROR, INPUT.IN_HY_ROR_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_HY_ROR.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    # basic_functions.create_assignment_matrix(INPUT.IN_NODE, INPUT.IN_HY_ROR, COLUMN_DEFINITIONS.COLUMNS_HY_ROR.node_id, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- Electric mobility A_E_MOB assignment matrix (Dimension: nodes x electirc mobilities)
    A_EMOB_RESIDENTIAL = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_EMOB_RESIDENTIAL, INPUT.IN_EMOB_RESIDENTIAL_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_EMOB.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id) #basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_E_MOB, INPUT.IN_E_MOB_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_E_MOB.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    A_EMOB_COMMERCIAL = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_EMOB_COMMERCIAL, INPUT.IN_EMOB_COMMERCIAL_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_EMOB.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    A_EMOB_QUICK = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_EMOB_QUICK, INPUT.IN_EMOB_QUICK_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_EMOB.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    #-- Heatpumps A_HP assignment matrix (Dimension: nodes x electirc mobilities)
    A_HP = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_NODE, INPUT.IN_HEATPUMP, INPUT.IN_HEATPUMP_ELEC_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_HEATPUMP.unit_name, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)

    #-- Power to Heat assignment matrix (Dimension: nodes x powerplants)
    A_PTH_E = basic_functions.create_assignment_matrix(INPUT.IN_NODE, INPUT.IN_PTH_ELEC, COLUMN_DEFINITIONS.COLUMNS_PTH.node_id, COLUMN_DEFINITIONS.COLUMNS_NODE.node_id)
    A_PTH_H = basic_functions.create_assignment_matrix(INPUT.IN_HEAT_NODE, INPUT.IN_PTH_HEAT, COLUMN_DEFINITIONS.COLUMNS_PTH.heatregion_id, COLUMN_DEFINITIONS.COLUMNS_HEATREGION.heatregion_id)
    #-- Heat only units assignment matrix (Dimension: nodes x powerplants)
    A_HEAT_ONLY = basic_functions.create_assignment_matrix(INPUT.IN_HEAT_NODE, INPUT.IN_HEAT_ONLY, COLUMN_DEFINITIONS.COLUMNS_HEAT_ONLY.heatregion_id, COLUMN_DEFINITIONS.COLUMNS_HEATREGION.heatregion_id)
    #-- Const heat load assignment matrix (Dimension: nodes x powerplants)
    A_CONST_HEAT_LOAD = basic_functions.create_assignment_matrix(INPUT.IN_HEAT_NODE, INPUT.IN_CONST_HEAT_LOAD, COLUMN_DEFINITIONS.COLUMNS_CONST_HEAT_LOAD.heatregion_id, COLUMN_DEFINITIONS.COLUMNS_HEATREGION.heatregion_id)
    #-- Heat demand assignment matrix (Dimension: nodes x powerplants)
    A_HEAT_DEMAND = basic_functions.create_assignment_matrix_timeseries(INPUT.IN_HEAT_NODE, INPUT.IN_HEAT_DEMAND, INPUT.IN_HEAT_DEMAND_TIMESERIES, COLUMN_DEFINITIONS.COLUMNS_HEAT_DEMAND.unit_name, COLUMN_DEFINITIONS.COLUMNS_HEATREGION.heatregion_id)
    # basic_functions.create_assignment_matrix(INPUT.IN_HEAT_NODE, INPUT.IN_HEAT_DEMAND, COLUMN_DEFINITIONS.COLUMNS_HEAT_DEMAND.heatregion_id, COLUMN_DEFINITIONS.COLUMNS_HEATREGION.heatregion_id)

    #--- Sets
    SET_PTG = 1:size(INPUT.IN_PTG,1)
    SET_GENS = 1:size(INPUT.IN_CONVENTIONAL,1)
    SET_RENS = 1:size(INPUT.IN_RENEWABLE,1)
    SET_LOADS = 1:size(INPUT.IN_LOAD,1)
    SET_BRANCHES = 1:size(INPUT.IN_BRANCHES,1)
    SET_NODES = 1:size(INPUT.IN_NODE,1)
    SET_TIME = 1:size(INPUT.IN_LOAD_TIMESERIES,1)
    SET_STORS = 1:size(INPUT.IN_STORAGE,1)
    SET_CHPS = 1:size(INPUT.IN_COMBINEDHEATPOWER_ELEC,1)
    SET_CHPS_BP = findall(x -> x=="CHP1", INPUT.IN_COMBINEDHEATPOWER_ELEC[!,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.operation_mode])
    SET_CHPS_EC = findall(x -> x=="CHP2", INPUT.IN_COMBINEDHEATPOWER_ELEC[!,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.operation_mode])
    SET_HEATREGIONS = 1:size(INPUT.IN_HEAT_NODE,1)
    SET_HVDC = 1:size(INPUT.IN_HVDC,1)
    SET_BIO = 1:size(INPUT.IN_BIO,1)
    SET_HY_ROR = 1:size(INPUT.IN_HY_ROR,1)
    SET_E_MOB_RESIDENTIAL = 1:size(INPUT.IN_EMOB_RESIDENTIAL,1)
    SET_E_MOB_COMMERCIAL = 1:size(INPUT.IN_EMOB_COMMERCIAL,1)
    SET_E_MOB_QUICK = 1:size(INPUT.IN_EMOB_QUICK,1)
    SET_HP = 1:size(INPUT.IN_HEATPUMP,1)
    SET_PTH = 1:size(INPUT.IN_PTH_ELEC,1)
    SET_HEAT_ONLY = 1:size(INPUT.IN_HEAT_ONLY,1)
    SET_CONST_HEAT_LOAD = 1:size(INPUT.IN_CONST_HEAT_LOAD,1)
    SET_HEAT_DEMAND = 1:size(INPUT.IN_HEAT_DEMAND,1)

#--- Create PTDF and matrices for voltage angle definition of power flow
    PTDF, DELTA_LF_MATRIX_LINE, DELTA_LF_MATRIX_NODE, A_FROM, A_TO, RESISTANCE = basic_functions.create_PTDF(INPUT.IN_BRANCHES, INPUT.IN_NODE, COLUMN_DEFINITIONS.COLUMNS_LINE,
        COLUMN_DEFINITIONS.COLUMNS_NODE, PARAMETER_SETTINGS.SLACK_NODE)

    P_MAX = INPUT.IN_BRANCHES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.Pmax]

    function get_loss_factor(RESISTANCE, P_MAX, rel_p, REF_VOLTAGE)
        LOSS_FACTOR_m = 2/3 .* RESISTANCE .* P_MAX .* (rel_p / REF_VOLTAGE^2)
        LOSS_FACTOR_b = -1/3 .* RESISTANCE .* P_MAX.^2 .* (rel_p^2 / REF_VOLTAGE^2)
        return LOSS_FACTOR_m, LOSS_FACTOR_b
    end

    # Loss factors
    LOSS_FACTOR_m, _ = get_loss_factor(RESISTANCE, P_MAX, 0.3125, PARAMETER_SETTINGS.REF_VOLTAGE)

    # LOSS_FACTOR_m0, LOSS_FACTOR_b0 = get_loss_factor(RESISTANCE, P_MAX, 0.675, PARAMETER_SETTINGS.REF_VOLTAGE)
    # LOSS_FACTOR_m1, LOSS_FACTOR_b1 = get_loss_factor(RESISTANCE, P_MAX, 0.1875, PARAMETER_SETTINGS.REF_VOLTAGE)
    # LOSS_FACTOR_m2, LOSS_FACTOR_b2 = get_loss_factor(RESISTANCE, P_MAX, 0.4, PARAMETER_SETTINGS.REF_VOLTAGE)
    # LOSS_FACTOR_m3, LOSS_FACTOR_b3 = get_loss_factor(RESISTANCE, P_MAX, 0.6375, PARAMETER_SETTINGS.REF_VOLTAGE)
    # LOSS_FACTOR_m4, LOSS_FACTOR_b4 = get_loss_factor(RESISTANCE, P_MAX, 0.8875, PARAMETER_SETTINGS.REF_VOLTAGE)
    # LOSS_FACTOR_m1a, LOSS_FACTOR_b1a = get_loss_factor(RESISTANCE, P_MAX, 0.18625, PARAMETER_SETTINGS.REF_VOLTAGE)
    # LOSS_FACTOR_m2a, LOSS_FACTOR_b2a = get_loss_factor(RESISTANCE, P_MAX, 0.7375, PARAMETER_SETTINGS.REF_VOLTAGE)
    # # Introducing an offset
    # LOSS_FACTOR_b2a = LOSS_FACTOR_b2a .- LOSS_FACTOR_b1a.*0.725
    # LOSS_FACTOR_b1a = LOSS_FACTOR_b1a .* 0

#--- Create strcuts for better function handling
    BRANCHES = (
        PTDF = PTDF,
        PMAX = INPUT.IN_BRANCHES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.Pmax],
        DELTA_LF_MATRIX_LINE = DELTA_LF_MATRIX_LINE,
        DELTA_LF_MATRIX_NODE = DELTA_LF_MATRIX_NODE,
        BRANCH_TYPE = INPUT.IN_BRANCHES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.Branch_type],
        VOLTAGE_F = INPUT.IN_BRANCHES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.Voltage_level_f],
        VOLTAGE_T = INPUT.IN_BRANCHES[:,COLUMN_DEFINITIONS.COLUMNS_LINE.Voltage_level_t],
        ASSIGNMENT_FROM = A_FROM,
        ASSIGNMENT_TO = A_TO,
        # RESISTANCE = RESISTANCE,
        # LOSS_FACTOR_m0 = LOSS_FACTOR_m0,
        # LOSS_FACTOR_b0 = LOSS_FACTOR_b0,
        # LOSS_FACTOR_m1 = LOSS_FACTOR_m1,
        # LOSS_FACTOR_b1 = LOSS_FACTOR_b1,
        # LOSS_FACTOR_m2 = LOSS_FACTOR_m2,
        # LOSS_FACTOR_b2 = LOSS_FACTOR_b2,
        # LOSS_FACTOR_m3 = LOSS_FACTOR_m3,
        # LOSS_FACTOR_b3 = LOSS_FACTOR_b3,
        # LOSS_FACTOR_m4 = LOSS_FACTOR_m4,
        # LOSS_FACTOR_b4 = LOSS_FACTOR_b4,
        # LOSS_FACTOR_m1a = LOSS_FACTOR_m1a,
        # LOSS_FACTOR_b1a = LOSS_FACTOR_b1a,
        # LOSS_FACTOR_m2a = LOSS_FACTOR_m2a,
        # LOSS_FACTOR_b2a = LOSS_FACTOR_b2a,
        LOSS_FACTOR_m = LOSS_FACTOR_m,
    )
    LOADS = (
        TIMESERIES = Matrix{Float64}(INPUT.IN_LOAD_TIMESERIES[:,2:end]),
    )
    RENEWABLES = (
        TIMESERIES = Matrix{Float64}(INPUT.IN_RENEWABLE_TIMESERIES[:,2:end]),
    )

    POWER_TO_GAS = (
        PMAX = INPUT.IN_PTG[:,COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.power_elec_max],
        ETA_P_MAX = INPUT.IN_PTG[:, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.eta_power_elec_max],
        ETA_P_MIN = INPUT.IN_PTG[:, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.eta_power_elec_min],
        H_MIN = INPUT.IN_PTG[:, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.heat_rate_elec_min],
        H_OP = INPUT.IN_PTG[:, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.heat_rate_marginal],
        R_G = INPUT.IN_PTG[:, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.min_load_factor],
        HYDROGEN_PRICE = INPUT.IN_PTG[:, COLUMN_DEFINITIONS.COLUMNS_POWERTOGAS.price_per_mwh],
    )
    CONVENTIONALS = (
        PMAX = INPUT.IN_CONVENTIONAL[:,COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.power_elec_max],
        FUEL_PRICE = INPUT.IN_CONVENTIONAL[:, COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.price_per_mwh],
        H_MIN = INPUT.IN_CONVENTIONAL[:, COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.heat_rate_elec_min],
        H_OP = INPUT.IN_CONVENTIONAL[:, COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.heat_rate_marginal],
        R_G = INPUT.IN_CONVENTIONAL[:, COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.min_load_factor],
        START_UP_COSTS = INPUT.IN_CONVENTIONAL[:, COLUMN_DEFINITIONS.COLUMNS_CONVENTIONAL.costs_startup],
    )
    STORAGES = (
        EMAX = INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.energy_elec_filling_level],
        EFF_CH = INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.eta_power_elec_charge],
        EFF_DIS = INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.eta_power_elec_discharge],
        SELF_DIS_RATIO = INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.self_dis_ratio],
        P_MAX_CH = INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.power_elec_max_charge],
        P_MAX_DIS = INPUT.IN_STORAGE[:,COLUMN_DEFINITIONS.COLUMNS_STORAGE.power_elec_max_discharge],
        PENALTY_LAST_TIMESTEP = PARAMETER_SETTINGS.PENALTY_LAST_TIMESTEP,
   )
    ASSIGNMENTS = (
        A_PTG = A_PTG,
        A_G = A_G,
        A_REN = A_REN,
        A_L  = A_L,
        A_STOR = A_STOR,
        A_CHP_EL = A_CHP_EL,
        A_CHP_H = A_CHP_H,
        A_HVDC = A_HVDC,
        A_HVDC_LOSS = A_HVDC_LOSS,
        A_BIO = A_BIO,
        A_HY_ROR = A_HY_ROR,
        A_EMOB_RESIDENTIAL = A_EMOB_RESIDENTIAL,
        A_EMOB_COMMERCIAL = A_EMOB_COMMERCIAL,
        A_EMOB_QUICK = A_EMOB_QUICK,
        A_HP = A_HP,
        A_PTH_E = A_PTH_E,
        A_PTH_H = A_PTH_H,
        A_HEAT_ONLY = A_HEAT_ONLY,
        A_CONST_HEAT_LOAD = A_CONST_HEAT_LOAD,
        A_HEAT_DEMAND = A_HEAT_DEMAND
    )
    SETS = (
        PTG = SET_PTG,
        GENS = SET_GENS,
        RENS = SET_RENS,
        LOADS = SET_LOADS,
        BRANCHES = SET_BRANCHES,
        NODES = SET_NODES,
        TIME = SET_TIME,
        STORS = SET_STORS,
        CHPS = SET_CHPS,
        CHPS_BP = SET_CHPS_BP,
        CHPS_EC = SET_CHPS_EC,
        HEATREGIONS = SET_HEATREGIONS,
        HVDCS = SET_HVDC,
        BIOS = SET_BIO,
        HY_ROR = SET_HY_ROR,
        EMOB_RESIDENTIAL = SET_E_MOB_RESIDENTIAL,
        EMOB_COMMERCIAL = SET_E_MOB_COMMERCIAL,
        EMOB_QUICK = SET_E_MOB_QUICK,
        HPS = SET_HP,
        PTH = SET_PTH,
        HEAT_ONLY = SET_HEAT_ONLY,
        CONST_HEAT_LOAD = SET_CONST_HEAT_LOAD,
        HEAT_DEMAND = SET_HEAT_DEMAND
    )
    PERIODS = (
        BASE_LENGTH = PARAMETER_SETTINGS.PERIOD_LENGTH,
        OVERLAP_LENGTH = PARAMETER_SETTINGS.PERIOD_OVERLAP_LENGTH,
        N_PERIODS = Int64(SET_TIME[end]/PARAMETER_SETTINGS.PERIOD_LENGTH) # NPeriods (number of periods)
    )
    COMBINEDHEATPOWERS = (
        PMAX = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.power_elec_max],
        QMAX = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.q_heat_max],
        FUEL_PRICE = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.price_per_mwh],
        H_MIN = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.heat_rate_elec_min],
        H_OP = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.heat_rate_marginal],
        R_G = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.min_load_factor],
        CB = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.chp_power_heat_ratio],
        CV = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.chp_power_loss],
        START_UP_COSTS = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.costs_startup],
        CHP_TYPE = INPUT.IN_COMBINEDHEATPOWER_ELEC[:,COLUMN_DEFINITIONS.COLUMNS_COMBINEDHEATPOWER.chp_mode]
    )
    HEAT_DEMAND = (
        TIMESERIES = Matrix{Float64}(INPUT.IN_HEAT_DEMAND_TIMESERIES[:,2:end]),
    )
    HVDCS = (
        P_MAX = INPUT.IN_HVDC[:,COLUMN_DEFINITIONS.COLUMNS_HVDC.Pmax],
        l0_const_loss = INPUT.IN_HVDC[:,COLUMN_DEFINITIONS.COLUMNS_HVDC.l0_const_loss],
        l1_proportional_loss = INPUT.IN_HVDC[:,COLUMN_DEFINITIONS.COLUMNS_HVDC.l1_proportional_loss]
    )
    BIOS = (
        PRICE = INPUT.IN_BIO[:,COLUMN_DEFINITIONS.COLUMNS_BIO.price_per_mwh],
        TIMESERIES = Matrix{Float64}(INPUT.IN_BIO_TIMESERIES[:,2:end])
    )
    HY_ROR = (
        TIMESERIES = Matrix{Float64}(INPUT.IN_HY_ROR_TIMESERIES[:,2:end]),
    )
    EMOB_RESIDENTIAL = (
        TIMESERIES = Matrix{Float64}(INPUT.IN_EMOB_RESIDENTIAL_TIMESERIES[:,2:end]),
        PMAX  = Matrix{Float64}(INPUT.IN_EMOB_RESIDENTIAL_TIMESERIES_PMAX[:,2:end])
    )

    if (PARAMETER_SETTINGS.FLEX_MARKET.EMOB == "Flex")
        @set EMOB_RESIDENTIAL.PMAX = max.(2*EMOB_RESIDENTIAL.TIMESERIES, EMOB_RESIDENTIAL.PMAX)
        println("EMOB load is limited to max. 2x grid friendly charging timeseries for the \"Flex\" model (system oriented smart charging).")
    end

    EMOB_COMMERCIAL = (
        TIMESERIES = Matrix{Float64}(INPUT.IN_EMOB_COMMERCIAL_TIMESERIES[:,2:end]),
    )
    EMOB_QUICK = (
        TIMESERIES = Matrix{Float64}(INPUT.IN_EMOB_QUICK_TIMESERIES[:,2:end]),
    )



    HEATPUMPS = (
        P_MAX_HP = INPUT.IN_HEATPUMP[:,COLUMN_DEFINITIONS.COLUMNS_HEATPUMP.power_elec_max],
        COP_TIMESERIES = Matrix{Float64}(INPUT.IN_HEATPUMP_COP_TIMESERIES[:,2:end]),
        P_MAX_HE = INPUT.IN_HEATPUMP[:, COLUMN_DEFINITIONS.COLUMNS_HEATPUMP.rod_power_elec_max],
        Q_MAX = INPUT.IN_HEATPUMP[:, COLUMN_DEFINITIONS.COLUMNS_HEATPUMP.energy_thermal_filling_level],
        NU_SELF_DIS = INPUT.IN_HEATPUMP[:, COLUMN_DEFINITIONS.COLUMNS_HEATPUMP.nu_thermal_self_discharge],
        ELEC_TIMESERIES = Matrix{Float64}(INPUT.IN_HEATPUMP_ELEC_TIMESERIES[:,2:end]),
        HEAT_TIMESERIES = Matrix{Float64}(INPUT.IN_HEATPUMP_HEAT_TIMESERIES[:,2:end]),
        PENALTY_LAST_TIMESTEP = PARAMETER_SETTINGS.PENALTY_LAST_TIMESTEP,
    )
    PTH = ( 
        PMAX = INPUT.IN_PTH_HEAT[:,COLUMN_DEFINITIONS.COLUMNS_PTH.power_elec_max],
    )
    HEAT_ONLY = (
        PMAX = INPUT.IN_HEAT_ONLY[:,COLUMN_DEFINITIONS.COLUMNS_HEAT_ONLY.power_heat_max],
        FUEL_PRICE = INPUT.IN_HEAT_ONLY[:,COLUMN_DEFINITIONS.COLUMNS_HEAT_ONLY.price_per_mwh]
    )
    CONST_HEAT_LOAD = (
        Q_MAX = INPUT.IN_CONST_HEAT_LOAD[:,COLUMN_DEFINITIONS.COLUMNS_CONST_HEAT_LOAD.power_heat_max],
    )

#--- preparing data for modelinput
    MODEL_INPUT = (
        BRANCHES = BRANCHES,
        LOADS = LOADS,
        RENEWABLES = RENEWABLES,
        CONVENTIONALS = CONVENTIONALS,
        STORAGES = STORAGES,
        ASSIGNMENTS = ASSIGNMENTS,
        SETS = SETS,
        PERIODS = PERIODS,
        COMBINEDHEATPOWERS = COMBINEDHEATPOWERS,
        HEAT_DEMAND = HEAT_DEMAND,
        HVDCS = HVDCS,
        BIOS = BIOS,
        HY_ROR = HY_ROR,
        EMOB_RESIDENTIAL = EMOB_RESIDENTIAL,
        EMOB_COMMERCIAL = EMOB_COMMERCIAL,
        EMOB_QUICK = EMOB_QUICK,
        HEATPUMPS = HEATPUMPS,
        PTH = PTH,
        HEAT_ONLY = HEAT_ONLY,
        CONST_HEAT_LOAD = CONST_HEAT_LOAD,
        PTG = POWER_TO_GAS
    )


    if PARAMETER_SETTINGS.TRANSFORMER_RESTRICTIONS == "NOT ACTIVE"
        # find considered branches
        idx_trafos = MODEL_INPUT.BRANCHES.BRANCH_TYPE .== "Transformer"
        print("Transformers: Amount = " * string(sum(idx_trafos)))
        # Select considered branches
        if sum(idx_trafos) > 0
            println(" --> For " * string(sum(idx_trafos)) * " transformers, thermal limit will be increased to: " * string(9999999) * ".")
            MODEL_INPUT.BRANCHES.PMAX[idx_trafos] = ones(sum(idx_trafos), 1) * 9999999
        end
    end
    
    if PARAMETER_SETTINGS.LINE_RESTRICTIONS_VOLTAGE > 0
        # find considered branches
        #print(MODEL_INPUT.BRANCHES.BRANCH_TYPE .== "Line" )
        idx_lines = MODEL_INPUT.BRANCHES.BRANCH_TYPE .== "Line"
        #println("lines: " * string(sum(idx_lines)))
        idx_low_voltage = ((MODEL_INPUT.BRANCHES.BRANCH_TYPE .== "Line") .& (MODEL_INPUT.BRANCHES.VOLTAGE_F .< PARAMETER_SETTINGS.LINE_RESTRICTIONS_VOLTAGE))
        print("Lines with voltage < " * string(PARAMETER_SETTINGS.LINE_RESTRICTIONS_VOLTAGE) * "kV: Amount = " * string(sum(idx_low_voltage)))
        # Select considered branches
        if sum(idx_low_voltage) > 0
            println(" --> For " * string(sum(idx_low_voltage)) * " lines, thermal limit will be increased to: " * string(999999999) * " to assume DSO grid is sufficiently dimensioned.")
            MODEL_INPUT.BRANCHES.PMAX[idx_low_voltage] = ones(sum(idx_low_voltage), 1) * 999999999
        end
    end

    if PARAMETER_SETTINGS.N_MINUS_1_FACTOR < 1
        # find considered branches
        #print(MODEL_INPUT.BRANCHES.BRANCH_TYPE .== "Line" )
        idx_lines = MODEL_INPUT.BRANCHES.BRANCH_TYPE .== "Line"
        #println("lines: " * string(sum(idx_lines)))
        #println("param_low_threshold" * string(PARAMETER_SETTINGS.N_MINUS_1_VOLTAGE))
        idx_high_voltage = ((MODEL_INPUT.BRANCHES.BRANCH_TYPE .== "Line") .& (MODEL_INPUT.BRANCHES.VOLTAGE_F .>= PARAMETER_SETTINGS.LINE_RESTRICTIONS_VOLTAGE))
        print("Lines with voltage >= " * string(PARAMETER_SETTINGS.LINE_RESTRICTIONS_VOLTAGE) * " kV: Amount = " * string(sum(idx_high_voltage)))
        # Select considered branches
        if sum(idx_high_voltage) > 0
            println(" --> For " * string(sum(idx_high_voltage)) * " lines, thermal limit will be multiplied with: " * string(PARAMETER_SETTINGS.N_MINUS_1_FACTOR) * " to depict \"n-1\" reserve.")
            MODEL_INPUT.BRANCHES.PMAX[idx_high_voltage] = PARAMETER_SETTINGS.N_MINUS_1_FACTOR .* MODEL_INPUT.BRANCHES.PMAX[idx_high_voltage]
        end
    end

    return MODEL_INPUT
end


function load_input_test_data(TEST_DATA_DEFINITIONS, PARAMETER_SETTINGS)
    #--- Load input
    FILE = TEST_DATA_DEFINITIONS.FILE
    SHEET = TEST_DATA_DEFINITIONS.SHEET
    CALC_TIME = PARAMETER_SETTINGS.CALC_TIME

    INPUT = (
    IN_NODE = DataFrame(XLSX.readtable(FILE, SHEET.node)...),
    IN_BRANCHES = DataFrame(XLSX.readtable(FILE, SHEET.line)...),
    IN_CONVENTIONAL = DataFrame(XLSX.readtable(FILE, SHEET.conventional)...),
    IN_RENEWABLE = DataFrame(XLSX.readtable(FILE, SHEET.renewable)...),
    IN_RENEWABLE_TIMESERIES = DataFrame(XLSX.readtable(FILE, SHEET.renewable_timeseries)...)[CALC_TIME,:],
    IN_LOAD = DataFrame(XLSX.readtable(FILE, SHEET.load)...),
    IN_LOAD_TIMESERIES = DataFrame(XLSX.readtable(FILE, SHEET.load_timeseries)...)[CALC_TIME,:],
    IN_STORAGE = DataFrame(XLSX.readtable(FILE, SHEET.storage)...),
    IN_COMBINEDHEATPOWER = DataFrame(XLSX.readtable(FILE, SHEET.chp)...),
    IN_HEATLOAD = DataFrame(XLSX.readtable(FILE, SHEET.heatload)...),
    IN_HEATLOAD_TIMESERIES = DataFrame(XLSX.readtable(FILE, SHEET.heatload_timeseries)...)[CALC_TIME,:],
    IN_HEATREGION = DataFrame(XLSX.readtable(FILE, SHEET.heatregion)...),
    IN_HVDC = DataFrame(XLSX.readtable(FILE, SHEET.HVDC)...),
    IN_BIO = DataFrame(XLSX.readtable(FILE, SHEET.bio)...),
    IN_BIO_TIMESERIES = DataFrame(XLSX.readtable(FILE, SHEET.bio_timeseries)...)[CALC_TIME,:],
    IN_E_MOB = DataFrame(XLSX.readtable(FILE, SHEET.e_mob)...),
    IN_E_MOB_TIMESERIES = DataFrame(XLSX.readtable(FILE, SHEET.e_mob_timeseries)...)[CALC_TIME,:],
    IN_HEATPUMP = DataFrame(XLSX.readtable(FILE, SHEET.heatpump)...),
    IN_HEATPUMP_TIMESERIES = DataFrame(XLSX.readtable(FILE, SHEET.heatpump_timeseries)...)[CALC_TIME,:],
    )
    return INPUT
end



function input_parameter_definitions_columns()
    #--- Column definitions
    COLUMNS = (
        COLUMNS_NODE = (
            node_id="node_id",
            node_name="name"
        ),
        COLUMNS_LINE = (
            name="name",
            node_id1="node_id1",
            node_id2="node_id2",
            voltage_level="voltage_level",
            Pmax="Pmax",
            x="x",
            length="length",
            imax="imax",
            Branch_type="Branch_type",
            Voltage_level_f= "Voltage_level_f",
            Voltage_level_t= "Voltage_level_t",
            r="r",
        ),
        COLUMNS_CONVENTIONAL = (
            unit_name="unit_name",
            node_id="node_id",
            power_elec_max="power_elec_max",
            power_elec_min="power_elec_min",
            eta_power_elec_max="eta_power_elec_max",
            eta_power_elec_min="eta_power_elec_min",
            heat_rate_elec_min="heat_rate_elec_min",
            heat_rate_marginal="heat_rate_marginal",
            min_load_factor="min_load_factor",
            price_per_mwh="price_per_mwh",
            costs_startup = "costs_startup",
        ),
        COLUMNS_POWERTOGAS = (
            unit_name="unit_name",
            node_id="node_id",
            power_elec_max="power_elec_max",
            power_elec_min="power_elec_min",
            eta_power_elec_max="eta_power_elec_max",
            eta_power_elec_min="eta_power_elec_min",
            heat_rate_elec_min="heat_rate_elec_min",
            heat_rate_marginal="heat_rate_marginal",
            min_load_factor="min_load_factor",
            price_per_mwh="price_per_mwh",
        ),
        COLUMNS_RENEWABLE = (
            unit_name="unit_name",
            node_id="node_id",
        ),
        COLUMNS_LOAD = (
            unit_name="unit_name",
            node_id="node_id",
        ),
        COLUMNS_STORAGE = (
            unit_name="unit_name",
            node_id="node_id",
            energy_elec_filling_level="energy_elec_filling_level",
            eta_power_elec_charge="eta_power_elec_charge",
            eta_power_elec_discharge="eta_power_elec_discharge",
            power_elec_max_discharge="power_elec_max_discharge",
            power_elec_max_charge="power_elec_max_charge",
            self_dis_ratio="self_dis_ratio",
        ),
        COLUMNS_COMBINEDHEATPOWER = (
            unit_name="unit_name",
            node_id="node_id",
            power_elec_max="power_elec_max",
            q_heat_max="q_heat_max",
            price_per_mwh="price_per_mwh",
            heat_rate_elec_min="heat_rate_elec_min",
            heat_rate_marginal="heat_rate_marginal",
            min_load_factor="min_load_factor",
            chp_power_heat_ratio="chp_power_heat_ratio",
            chp_power_loss="chp_power_loss",
            operation_mode="operation_mode",
            heatregion_id="heat_node_id",
            power_elec_min="power_elec_min",
            eta_power_elec_min="eta_power_elec_min",
            eta_power_elec_max="eta_power_elec_max",
            costs_startup="costs_startup",
            chp_mode="operation_mode",
        ),
        COLUMNS_HEAT_DEMAND = (
            unit_name="unit_name",
            heatregion_id="heat_node_id",
        ),
        COLUMNS_PTH = (
            unit_name="unit_name",
            node_id = "node_id",
            heatregion_id="heat_node_id",
            power_elec_max = "power_elec_max",
        ),
        COLUMNS_HEAT_ONLY = (
            unit_name="unit_name",
            heatregion_id="heat_node_id",
            power_heat_max = "power_heat_max",
            price_per_mwh = "price_per_mwh",
        ),
        COLUMNS_CONST_HEAT_LOAD = (
            unit_name="unit_name",
            heatregion_id="heat_node_id",
            power_heat_max = "power_heat_max",
        ),
        COLUMNS_HEATREGION = (
            heatregion_id="heat_node_id",
            power_heat_max = "power_heat_max",
        ),
        COLUMNS_HVDC = (
            unit_name="name",
            node_id1="node_id1",
            node_id2="node_id2",
            l0_const_loss="l0_const_loss",
            l1_proportional_loss="l1_proportional_loss",
            Pmax="maxcap",
        ),
        COLUMNS_BIO = (
            unit_name="unit_name",
            node_id="node_id",
            price_per_mwh="price_per_mwh",
        ),
        COLUMNS_HY_ROR = (
            unit_name="unit_name",
            node_id="node_id",
        ),
        COLUMNS_EMOB = (
            unit_name="unit_name",
            node_id="node_id",
        ),
        COLUMNS_HEATPUMP = (
            unit_name="unit_name",
            node_id="node_id",
            power_elec_max = "heat_pump_power_elec_max",
            rod_power_elec_max = "heat_pump_rod_power_elec_max",
            energy_thermal_filling_level = "heat_pump_energy_thermal_filling_level",
            nu_elec_self_discharge = "heat_pump_nu_elec_self_discharge",
            energy_elec_filling_level = "heat_pump_energy_elec_filling_level",
            eta_power_elec_charge = "heat_pump_eta_power_elec_charge",
            eta_power_elec_discharge = "heat_pump_eta_power_elec_discharge",
            nu_thermal_self_discharge = "heat_pump_nu_thermal_self_discharge",
        ),

        COLUMNS_TRANSFORMER = (
            unit_name="unit_name",
            node_id1="node_id1",
            node_id2="node_id2",
            short_circuit_voltage="short_circuit_voltage",
            nominal_apparent_power="nominal_apparent_power",
        ),

        COLUMNS_BEV = (
            unit_name="unit_name",
            node_id="node_id"
        )
    )
    return COLUMNS
end

end
