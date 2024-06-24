module load_data

using CSV, DataFrames


function load_file_paths()
    input_data_folder_path = rsplit(dirname(@__FILE__),"\\"; limit=2)[1] * "\\input_data\\"

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
        return_tuple = merge(return_tuple, [file_list[ct]=>string(input_data_folder_path * FILE_PATHS[ct,2]),])

    end

    return return_tuple

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

    if (PARAMETER_SETTINGS.FLEX_MARKET.EMOB == "Flex")
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
