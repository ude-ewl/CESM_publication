module setup_functions

using Pkg

function add_additional_packages(str_additional_packages)

    # depending on your Gurobi and CPLEX standalone solver version you may need to up-/downgrade julia's Gurobi and CPLEX library, e.g.
    # Pkg.add(Pkg.PackageSpec(;name="CPLEX", version="0.7.8"))

    if occursin("CPLEX", str_additional_packages)
        Pkg.add("CPLEX")
        Pkg.build("CPLEX")
    end
    if occursin("Gurobi", str_additional_packages)
        Pkg.add("Gurobi")
        Pkg.build("Gurobi")
    end
end

function create_folders()

    # check if data folders exist, else indicate what to do
    errorsum = 0
    folderstring = rsplit(dirname(@__FILE__),"\\"; limit=2)[1] * "\\input_data"

    if ~isdir(folderstring)
        mkdir(folderstring)
        errorsum = errorsum + 1
    end

    if errorsum > 0
        @assert(false, "Couldn't locate the folder " * folderstring * ". Execution is terminated as input data is missing. An empty folder \"input_data\" was created. Please fill this folder with the input data.") 
    end
end

end
