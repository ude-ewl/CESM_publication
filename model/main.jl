# activate the project
cd(dirname(@__FILE__))
using Pkg

Pkg.activate("model_plain_vanilla")

include("setup_functions.jl")
if length(Pkg.installed()) == 0 # if no packages are in environment then install all packages as in Project.toml
    # download all the libraries as in Manifest.toml
    Pkg.instantiate()
    add_additional_packages = "CPLEX"
    setup_functions.add_additional_packages(add_additional_packages)
end 

setup_functions.create_folders()


# import the necessary packages
using JuMP, XLSX, CSV, DataFrames, Setfield, LinearAlgebra, Dates, StatsBase


# ======
#-- Load choice file
# ======
include("choice_file.jl")
PARAMETER_SETTINGS = choice_file.input_parameter_settings();

include("model_chains.jl")
model_chains.calling_solver_package(PARAMETER_SETTINGS)


# ======
# ==== load data
# ======
include("input_functions.jl")
#-- load column definitions (names) for tables
COLUMN_DEFINITIONS = input_functions.input_parameter_definitions_columns();
#-- Load input data from csv files and preprocess
INPUT_RAW = input_functions.load_input_data(COLUMN_DEFINITIONS, PARAMETER_SETTINGS);
#--- Restructure model input (MI) data for model input
INPUT = input_functions.fill_technology_structs(INPUT_RAW, COLUMN_DEFINITIONS, PARAMETER_SETTINGS);


# ======
# ==== arbitrary input data manipulations to test the model
# ======
#INPUT = @set INPUT.RENEWABLES.TIMESERIES =  0.001 * INPUT.RENEWABLES.TIMESERIES; # Setfield workaround, as struct SET is principally immuable (for increased performance)
#INPUT = @set INPUT.COMBINEDHEATPOWERS.PMAX =  0.001 * INPUT.COMBINEDHEATPOWERS.PMAX; # Setfield workaround, as struct SET is principally immuable (for increased performance)
#INPUT = @set INPUT.CONVENTIONALS.FUEL_PRICE =  LinRange(50, 100, 205)
INPUT = @set INPUT.BRANCHES.PMAX =  0.5 * INPUT.BRANCHES.PMAX


# ======
# ==== run the model
# ======

RESULTS_ALL = model_chains.rolling_optimization(INPUT, PARAMETER_SETTINGS)


# ======
# ==== save the results
# ======
what_to_save = Dict("aggregated_xlsx" => true, "detailed_csv" => false, "online_model" => false, "AoTRF" => false, "jld2" => false)
include("output_save_functions.jl")
output_save_functions.save_output_combination(RESULTS_ALL, INPUT, INPUT_RAW, COLUMN_DEFINITIONS, what_to_save, PARAMETER_SETTINGS)


# ======
# ==== plotting or results
# ======
what_to_plot = Dict("date" => "2023_02_11_15_55_44" , "rd_date" => "2023_02_11_15_55_45", "LMP_aggregated" => false, "nodal_prices" => false, "line_loading"=>false, "rd_quantities" => true)
include("output_plot_functions.jl")
output_plot_functions.plot_selection(what_to_plot, INPUT, RESULTS_ALL, RD_RESULTS_ALL, PARAMETER_SETTINGS)