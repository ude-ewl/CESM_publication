module load_data

using CSV, DataFrames






function load_file_paths()
    tmp_data_folder_path = rsplit(dirname(@__FILE__),"\\"; limit=2)[1] * "\\tmp_data\\"

    FILE_PATHS = CSV.read(raw"FILE_PATHS.csv",DataFrame)

    file_list = [:local_path, :elec_ptg_hydrogen, :assignment_elec_ptg_hydrogen, :elec_demand, :elec_pv, :elec_wind,
    :elec_bio, :elec_hy_ror, :elec_bev_residential, :elec_bev_commercial_combined, :elec_bev_quickcharge,
    :elec_hp_cop, :elec_hp_heat_demand, :assignment_elec_demand, :assignment_elec_pv, :assignment_elec_wind,
    :assignment_elec_bio, :assignment_elec_hy_ror, :assignment_elec_bev_residential, :assignment_elec_bev_commercial_combined, :assignment_elec_bev_quickcharge,
    :elec_hp_parameters, :assignment_elec_hp, :chp_parameters, :assignment_chp_units_elec, :assignment_chp_units_heat, 
    :el_parameters, :assignment_el_units, :el_storage_parameters, :assignment_el_storage_units, :elec_network_nodes,
    :elec_network_ac_lines, :elec_network_hvdc_lines, :elec_network_transformers, :heat_network_nodes, :heat_pth,
    :heat_only, :const_heat_load, :assignment_elec_pth, :assignment_heat_pth, :assignment_heat_only,
    :assignment_const_heat_load, :heat_demand, :assignment_heat_demand, :assignment_elec_bev_residential_pmax
    ]

    return_tuple = NamedTuple()

    for ct = collect(1:size(file_list)[1])
        return_tuple = merge(return_tuple, [file_list[ct]=>string(tmp_data_folder_path * FILE_PATHS[ct,2]),])

    end

    return return_tuple

end

function save_data_from_db(DB_CONNECTION)
    ### +++++++++++++++++++++++++++ Update local files +++++++++++++++++++++++++++
    
    ## Power-to-Gas
    # Get file_path and local_path of Power-to-gas Hydrogen paramters
    file_path_elec_ptg_hydrogen_parameters, local_path = database_functions.download_ptg_hydrogen_parameters(DB_CONNECTION)
    # Get file_path and local_path ofPower-to-gas Hydrogen assignment data
    file_path_assignment_elec_ptg_hydrogen, local_path = database_functions.download_unit_2_elec_node_ptg_hydrogen_units(DB_CONNECTION)

    ## Electric vehicles
    # Get file_path and local_path of electricty demand
    file_path_bev_commercial_combined, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "EV", "BEV_commercial_combined")
    file_path_bev_quickcharge, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "EV", "BEV_quickcharge")
    file_path_bev_residential, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "EV", "BEV_residential")
    # Get file_path and local_path of electricty demand assignment data
    file_path_assignment_bev_commercial_combined, local_path = database_functions.download_unit_2_elec_node_elec_bev_commercial_combined(DB_CONNECTION)
    file_path_assignment_bev_quickcharge, local_path = database_functions.download_unit_2_elec_node_elec_bev_quickcharge(DB_CONNECTION)
    file_path_assignment_bev_residential, local_path = database_functions.download_unit_2_elec_node_elec_bev_residential(DB_CONNECTION)

    ## Electricty network
    # Get file_path and local_path of electricty network components
    file_path_elec_network_nodes, local_path = database_functions.download_elec_network_nodes(DB_CONNECTION)
    file_path_elec_network_ac_lines, local_path = database_functions.download_elec_network_ac_lines(DB_CONNECTION)
    file_path_elec_network_hvdc_lines, local_path = database_functions.download_elec_network_hvdc_lines(DB_CONNECTION)
    file_path_elec_network_transformers, local_path = database_functions.download_elec_network_transformers(DB_CONNECTION)

    ## Heat regions
    # Get heat nodes
    file_path_heat_network_nodes, local_path = database_functions.download_heat_network_nodes(DB_CONNECTION)
    # Get file_path and local_path of chp unit paramters
    file_path_pth_parameters, local_path = database_functions.download_pth_parameters(DB_CONNECTION)
    file_path_heat_only_parameters, local_path = database_functions.download_heat_only_parameters(DB_CONNECTION)
    file_path_const_heat_load_parameters, local_path = database_functions.download_const_heat_load_parameters(DB_CONNECTION)
    # Get file_path and local_path of chp unit assignment data
    file_path_assignment_elec_pth_units, local_path = database_functions.download_unit_2_elec_node_pth_units(DB_CONNECTION)
    file_path_assignment_heat_pth_units, local_path = database_functions.download_unit_2_heat_node_pth_units(DB_CONNECTION)
    file_path_assignment_heat_only_units, local_path = database_functions.download_unit_2_heat_node_heat_only_units(DB_CONNECTION)
    file_path_assignment_const_heat_load_units, local_path = database_functions.download_unit_2_heat_node_const_heat_load_units(DB_CONNECTION)

    ## Heat demand
    # Get file_path and local_path of electricty demand
    file_path_heat_demand, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "HEAT_DEM", "heat_demand")
    # Get file_path and local_path of electricty demand assignment data
    file_path_assignment_heat_demand, local_path = database_functions.download_unit_2_heat_node_heat_demand(DB_CONNECTION)

    ## CHP units
    # Get file_path and local_path of chp unit paramters
    file_path_chp_parameters, local_path = database_functions.download_chp_parameters(DB_CONNECTION)
    # Get file_path and local_path of data for chp unit to elec node assignments
    file_path_assignment_chp_units_elec, local_path = database_functions.download_unit_2_elec_node_chp_units(DB_CONNECTION)
    # Get file_path and local_path of data for chp unit to heat node assignments
    file_path_assignment_chp_units_heat, local_path = database_functions.download_unit_2_heat_node_chp_units(DB_CONNECTION)

    ## Electricty demand
    # Get file_path and local_path of electricty demand
    file_path_elec_demand, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "ELEC_DEM", "elec_load")
    # Get file_path and local_path of electricty demand assignment data
    file_path_assignment_elec_demand, local_path = database_functions.download_unit_2_elec_node_elec_demand(DB_CONNECTION)

    ## Renewables
    # Get file_path and local_path of renewable electricity generation
    file_path_elec_pv, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "PV", "PV")
    file_path_elec_wind, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "Wind_On", "WIND")
    file_path_elec_bio, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "BIO", "BIO")
    file_path_elec_hy_ror, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "HY_ROR", "HY_ROR")
    ## EMOB
    file_path_elec_emob_res, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "EV", "bev_residential")
    file_path_elec_emob_comm, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "EV", "bev_commercial_combined")
    file_path_elec_emob_quick, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "EV", "bev_quickcharge")
  
    # Get file_path and local_path of renewable electricity generation data
    file_path_assignment_elec_pv, local_path = database_functions.download_unit_2_elec_node_elec_pv(DB_CONNECTION)
    file_path_assignment_elec_wind, local_path = database_functions.download_unit_2_elec_node_elec_wind(DB_CONNECTION)
    file_path_assignment_elec_bio, local_path = database_functions.download_unit_2_elec_node_elec_bio(DB_CONNECTION)
    file_path_assignment_elec_hy_ror, local_path = database_functions.download_unit_2_elec_node_elec_hy_ror(DB_CONNECTION)
    
    file_path_assignment_elec_bev_residential, local_path = database_functions.download_unit_2_elec_node_elec_bev_residential(DB_CONNECTION)
    file_path_assignment_elec_bev_commercial, local_path = database_functions.download_unit_2_elec_node_elec_bev_commercial_combined(DB_CONNECTION)
    file_path_assignment_elec_bev_quick, local_path = database_functions.download_unit_2_elec_node_elec_bev_quickcharge(DB_CONNECTION)

    ## Heat pumps
    # Get file_path and local_path of heat pump's COP and heat demand
    file_path_elec_hp_cop, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "HP_ELEC_DEM", "hp_cop")
    file_path_elec_hp_heat_demand, local_path = database_functions.download_timeseries_data(DB_CONNECTION, "HP_ELEC_DEM", "hp_demand")
    # Get file_path and local_path of heat pump unit paramters
    file_path_elec_hp_parameters, local_path = database_functions.download_el_heat_pump_parameters(DB_CONNECTION)
    # Get file_path and local_path of heat pump's unit assignment data
    file_path_assignment_elec_hp, local_path = database_functions.download_unit_2_elec_node_elec_heat_pump(DB_CONNECTION)

    ## Electricty only units
    # Get file_path and local_path of electricty only unit paramters
    file_path_el_parameters, local_path = database_functions.download_el_parameters(DB_CONNECTION)
    # Get file_path and local_path of electricty only unit assignment data
    file_path_assignment_el_units, local_path = database_functions.download_unit_2_elec_node_elec_only_units(DB_CONNECTION)

    ## Storage units
    # Get file_path and local_path of electricty only unit paramters
    file_path_el_storage_parameters, local_path = database_functions.download_el_storage_parameters(DB_CONNECTION)
    # Get file_path and local_path of electricty only unit assignment data
    file_path_assignment_el_storage_units, local_path = database_functions.download_unit_2_elec_node_el_storage_units(DB_CONNECTION)
    FILE_PATH = Dict("file_path_elec_ptg_hydrogen_parameters", file_path_elec_ptg_hydrogen_parameters)

    ## Aggregate results in one struct
    FILE_PATH = structs.f_file_paths(local_path,
    	# Power-to-gas Hydrogen parameters and assignment
        file_path_elec_ptg_hydrogen_parameters, file_path_assignment_elec_ptg_hydrogen,
        # Timeseries data of electricty demand, PV infeed, Wind infeed and biomass/biogas infeed
        file_path_elec_demand, file_path_elec_pv, file_path_elec_wind, file_path_elec_bio, file_path_elec_hy_ror,
        # Timeseries emob
        file_path_elec_emob_res, file_path_elec_emob_comm, file_path_elec_emob_quick,
        # Timeseries data of heat pump's COP and heat demand
        file_path_elec_hp_cop, file_path_elec_hp_heat_demand,
        # Assignment data of electricty demand, PV infeed, wind infeed, biomass/biogas infeed and hydro run-of-river (hy_ror)
        file_path_assignment_elec_demand, file_path_assignment_elec_pv, file_path_assignment_elec_wind,
        file_path_assignment_elec_bio, file_path_assignment_elec_hy_ror,
        file_path_assignment_elec_bev_residential, file_path_assignment_elec_bev_commercial, file_path_assignment_elec_bev_quick,
        # Parameters heat pump units
        file_path_elec_hp_parameters, file_path_assignment_elec_hp,
        # Parameters and assignment data of chp units
        file_path_chp_parameters, file_path_assignment_chp_units_elec, file_path_assignment_chp_units_heat,
        # Parameters and assignment data of electricty-only units
        file_path_el_parameters, file_path_assignment_el_units,
        # Parameters and assignment data of electricity storages
        file_path_el_storage_parameters, file_path_assignment_el_storage_units,
        # Network data
        file_path_elec_network_nodes, file_path_elec_network_ac_lines, file_path_elec_network_hvdc_lines, file_path_elec_network_transformers,
        # Heat region nodes
        file_path_heat_network_nodes,
        # Heat region unit parameters
        file_path_pth_parameters, file_path_heat_only_parameters, file_path_const_heat_load_parameters,
        # Heat region assignment data
        file_path_assignment_elec_pth_units, file_path_assignment_heat_pth_units, file_path_assignment_heat_only_units, file_path_assignment_const_heat_load_units,
        # Heat demand timeseries data and assignments
        file_path_heat_demand, file_path_assignment_heat_demand
        )

    FILE_PATHS = DataFrame(Name=collect(fieldnames(typeof(FILE_PATH))), Path=[n=getfield(FILE_PATH, n) for n in fieldnames(typeof(FILE_PATH))])
    CSV.write("FILE_PATHS.csv", FILE_PATHS)

    return FILE_PATH
end

function load_data_from_tmp_folder(PARAMETER_SETTINGS)

    FILE_PATH = load_file_paths()
    
    ### +++++++++++++++++++++++++++ Load local files +++++++++++++++++++++++++++
    ## Power-to-Gas
    # Load Power-to-gas Hydrogen parameters as DataFrame
    data_elec_ptg_hydrogen = DataFrame(CSV.read(FILE_PATH.elec_ptg_hydrogen, DataFrame))
    # Load Power-to-gas Hydrogen assignment data as DataFrame
    assignment_elec_ptg_hydrogen = DataFrame(CSV.read(FILE_PATH.assignment_elec_ptg_hydrogen, DataFrame))

    ## Electricty demand
    # Load data for electricty demand as DataFrame
    data_elec_demand = DataFrame(CSV.read(FILE_PATH.elec_demand, DataFrame))
    # Load assignment data for electricty demand as DataFrame
    assignment_elec_demand = DataFrame(CSV.read(FILE_PATH.assignment_elec_demand, DataFrame))

    ## Renewables
    # Load data for electricty generation time-series as DataFrame
    data_elec_pv = DataFrame(CSV.read(FILE_PATH.elec_pv, DataFrame))
    data_elec_wind = DataFrame(CSV.read(FILE_PATH.elec_wind, DataFrame))
    data_elec_bio = DataFrame(CSV.read(FILE_PATH.elec_bio, DataFrame))
    data_elec_hy_ror = DataFrame(CSV.read(FILE_PATH.elec_hy_ror, DataFrame))
   
    # Emob unit2node
    data_elec_bev_residential = DataFrame(CSV.read(FILE_PATH.elec_bev_residential, DataFrame))
    data_elec_bev_commercial = DataFrame(CSV.read(FILE_PATH.elec_bev_commercial_combined, DataFrame))
    data_elec_bev_quick = DataFrame(CSV.read(FILE_PATH.elec_bev_quickcharge, DataFrame))

    if (PARAMETER_SETTINGS.FLEX_MARKET.EMOB == "Flex") | (PARAMETER_SETTINGS.FLEX_RD.EMOB == "Flex")
        data_elec_bev_residential_pmax = DataFrame(CSV.read(FILE_PATH.assignment_elec_bev_residential_pmax, DataFrame))
    else
        data_elec_bev_residential_pmax = data_elec_bev_residential
    end

    # Load assignment data for electricty generation as DataFrame
    assignment_elec_pv = DataFrame(CSV.read(FILE_PATH.assignment_elec_pv, DataFrame))
    assignment_elec_wind = DataFrame(CSV.read(FILE_PATH.assignment_elec_wind, DataFrame))
    assignment_elec_bio = DataFrame(CSV.read(FILE_PATH.assignment_elec_bio, DataFrame))
    assignment_elec_hy_ror = DataFrame(CSV.read(FILE_PATH.assignment_elec_hy_ror, DataFrame))
    assignment_elec_bev_residential = DataFrame(CSV.read(FILE_PATH.assignment_elec_bev_residential, DataFrame))
    assignment_elec_bev_commercial = DataFrame(CSV.read(FILE_PATH.assignment_elec_bev_commercial_combined, DataFrame))
    assignment_elec_bev_quick = DataFrame(CSV.read(FILE_PATH.assignment_elec_bev_quickcharge, DataFrame))

    ## Heat pumps
    # Load data of heat pump's COP and heat demand as DataFrame
    data_elec_hp_cop = DataFrame(CSV.read(FILE_PATH.elec_hp_cop, DataFrame))
    data_elec_hp_heat_demand = DataFrame(CSV.read(FILE_PATH.elec_hp_heat_demand, DataFrame))
    # Load heat pump's unit paramters as DataFrame
    data_elec_hp_parameters = DataFrame(CSV.read(FILE_PATH.elec_hp_parameters, DataFrame))
    # Load assignment data for electricty demand of heat pumps as DataFrame
    assignment_elec_hp = DataFrame(CSV.read(FILE_PATH.assignment_elec_hp, DataFrame))

    ## CHP units
    # Load chp unit paramters as DataFrame
    data_chp_parameters = DataFrame(CSV.read(FILE_PATH.chp_parameters, DataFrame))
    # Load assignment data for chp unit' to elec_nodes as DataFrame
    assignment_chp_units_elec = DataFrame(CSV.read(FILE_PATH.assignment_chp_units_elec, DataFrame))
    # Load assignment data for chp unit' to heat_nodes as DataFrame
    assignment_chp_units_heat = DataFrame(CSV.read(FILE_PATH.assignment_chp_units_heat, DataFrame))

    ## Electricty only units
    # Load electricty only unit parameters as DataFrame
    data_el_parameters = DataFrame(CSV.read(FILE_PATH.el_parameters, DataFrame))
    # Load electricty only unit parameters as DataFrame
    assignment_el_units  = DataFrame(CSV.read(FILE_PATH.assignment_el_units, DataFrame))

    ## Electricty only units
    # Load electricty only unit parameters as DataFrame
    data_el_storage_parameters = DataFrame(CSV.read(FILE_PATH.el_storage_parameters, DataFrame))
    # Load electricty only unit parameters as DataFrame
    assignment_el_storage_units  = DataFrame(CSV.read(FILE_PATH.assignment_el_storage_units, DataFrame))

    ## Electricty network
    # Load electricty network components as DataFrame
    data_elec_network_nodes = DataFrame(CSV.read(FILE_PATH.elec_network_nodes, DataFrame))
    data_elec_network_ac_lines = DataFrame(CSV.read(FILE_PATH.elec_network_ac_lines, DataFrame))
    data_elec_network_hvdc_lines = DataFrame(CSV.read(FILE_PATH.elec_network_hvdc_lines, DataFrame))
    data_elec_network_transformers = DataFrame(CSV.read(FILE_PATH.elec_network_transformers, DataFrame))

    ## Heat regions
    # Load heat region data
    data_heat_network_nodes = DataFrame(CSV.read(FILE_PATH.heat_network_nodes, DataFrame))
    # Load unit parameters as DataFrames
    data_heat_pth = DataFrame(CSV.read(FILE_PATH.heat_pth, DataFrame))
    data_heat_only = DataFrame(CSV.read(FILE_PATH.heat_only, DataFrame))
    data_const_heat_load = DataFrame(CSV.read(FILE_PATH.const_heat_load, DataFrame))
    data_heat_demand = DataFrame(CSV.read(FILE_PATH.heat_demand, DataFrame))
    # Load assignment data for heat region units as DataFrame
    assignment_elec_pth = DataFrame(CSV.read(FILE_PATH.assignment_elec_pth, DataFrame))
    assignment_heat_pth = DataFrame(CSV.read(FILE_PATH.assignment_heat_pth, DataFrame))
    assignment_heat_only = DataFrame(CSV.read(FILE_PATH.assignment_heat_only, DataFrame))
    assignment_const_heat_load = DataFrame(CSV.read(FILE_PATH.assignment_const_heat_load, DataFrame))
    assignment_heat_demand = DataFrame(CSV.read(FILE_PATH.assignment_heat_demand, DataFrame))

    data_timeseries = (elec_demand = data_elec_demand, 
                        elec_pv = data_elec_pv, 
                        elec_wind = data_elec_wind, 
                        elec_bio = data_elec_bio,
                        elec_hy_ror = data_elec_hy_ror,
                        elec_bev_residential = data_elec_bev_residential, 
                        elec_bev_commercial = data_elec_bev_commercial, 
                        elec_bev_quick = data_elec_bev_quick, 
                        elec_hp_cop = data_elec_hp_cop, 
                        elec_hp_heat_demand = data_elec_hp_heat_demand, 
                        heat_demand = data_heat_demand,
                        elec_bev_residential_pmax = data_elec_bev_residential_pmax
    )


    data_unit_parameters = (
        ptg_hydrogen_units = data_elec_ptg_hydrogen, 
        chp_units = data_chp_parameters, 
        elec_only_units = data_el_parameters, 
        el_storage_units = data_el_storage_parameters, 
        elec_hp = data_elec_hp_parameters,
        heat_pth = data_heat_pth, 
        heat_only = data_heat_only, 
        const_heat_load = data_const_heat_load
        )

    data_assignment = (
        ptg_hydrogen_units = assignment_elec_ptg_hydrogen, 
        elec_demand = assignment_elec_demand, 
        elec_pv = assignment_elec_pv,
        elec_wind = assignment_elec_wind, 
        elec_bio = assignment_elec_bio, 
        elec_hy_ror = assignment_elec_hy_ror, 
        elec_hp = assignment_elec_hp,
        elec_bev_residential = assignment_elec_bev_residential, 
        elec_bev_commercial = assignment_elec_bev_commercial, 
        elec_bev_quick = assignment_elec_bev_quick,
        chp_units_elec = assignment_chp_units_elec, 
        chp_units_heat = assignment_chp_units_heat, 
        elec_only_units = assignment_el_units, 
        el_storage_units = assignment_el_storage_units,
        elec_pth = assignment_elec_pth, 
        heat_pth = assignment_heat_pth, 
        heat_only = assignment_heat_only, 
        const_heat_load = assignment_const_heat_load, 
        heat_demand = assignment_heat_demand
        )

    data_elec_network = (
        nodes = data_elec_network_nodes, 
        ac_lines = data_elec_network_ac_lines, 
        hvdc_lines = data_elec_network_hvdc_lines, 
        transformers = data_elec_network_transformers
    )

    data_heat_network = (
        nodes = data_heat_network_nodes, 
    )

    return data_timeseries, data_unit_parameters, data_assignment, data_elec_network, data_heat_network
end


end
