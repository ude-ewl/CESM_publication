# ======
# === Setting up the environment
# ======
cd(dirname(@__FILE__))
using Pkg

# activate the julia environment given the informaiton in the Project.toml
Pkg.activate(".") # this folder

my_solver = "Gurobi" # "CPLEX"

include("setup_functions.jl")
if length(Pkg.installed()) == 0 
    Pkg.instantiate() # download all the packages
    setup_functions.add_additional_packages(my_solver) # Solvers must not only be installed, but also be built
    Pkg.build("Plots")
end 
setup_functions.create_folders() # sanity check if "input_data" folder exists

# import the necessary packages for the main.jl script
using Setfield


# ======
# === Load choice file
# ======
include("choice_file.jl")
PARAMETER_SETTINGS = choice_file.input_parameter_settings();

include("model_chains.jl")
model_chains.calling_solver_package(PARAMETER_SETTINGS)


# ======
# === load data
# ======
include("input_functions.jl")
#-- load column definitions (names) for tables
COLUMN_DEFINITIONS = input_functions.input_parameter_definitions_columns();
#-- Load input data from csv files and preprocess
INPUT_RAW = input_functions.load_input_data(COLUMN_DEFINITIONS, PARAMETER_SETTINGS);
#--- Restructure model input (MI) data for model input
INPUT = input_functions.fill_technology_structs(INPUT_RAW, COLUMN_DEFINITIONS, PARAMETER_SETTINGS);

bahn.de/aktuell

# ======
# ==== arbitrary input data manipulations to test the model
# ======
#INPUT = @set INPUT.RENEWABLES.TIMESERIES =  0.001 * INPUT.RENEWABLES.TIMESERIES;
#INPUT = @set INPUT.COMBINEDHEATPOWERS.PMAX =  0.001 * INPUT.COMBINEDHEATPOWERS.PMAX;
#INPUT = @set INPUT.CONVENTIONALS.FUEL_PRICE =  LinRange(50, 100, 205)
#INPUT = @set INPUT.BRANCHES.PMAX =  0.5 * INPUT.BRANCHES.PMAX


# ======
# ==== run the model
# ======

RESULTS_ALL = model_chains.rolling_optimization(INPUT, PARAMETER_SETTINGS)


# ======
# ==== save the results
# ======
include("output_save_functions.jl")
output_save_functions.save_aggregated_results_as_excel(RESULTS_ALL, INPUT, INPUT_RAW, PARAMETER_SETTINGS)

# ======
# ==== plotting or results
# ======
include("output_plot_functions.jl")
output_plot_functions.plot_lmp_aggregated(INPUT, RESULTS_ALL)
