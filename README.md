[![DOI](https://zenodo.org/badge/687057711.svg)](https://zenodo.org/doi/10.5281/zenodo.11396393)

# Cellular Energy System Model 
This repo contains the Cellular Energy System Model (CESM). The model can be used to perform reserach on market design of decarbonized and decrentralized energy systems.
The CESM model was established for the publication: 
```Schinke-Nendza, Aiko; Kramer, Hendrik; Khalid, Abuzar; Flatter, Felix; Rasti, Sasan; Uhlemeyer, Björn; Weber, Christoph: Modeling of large decentral multi-energy systems - An application to flexible heat pumps in a decarbonized scenario, TechRxiv, 2024. DOI:10.36227/techrxiv.170861965.51606148/v1```


# Repo Structure
* ```input_data```: This is the folder of the input data. Data is available at the following data reposity hosted by Zenodo (see below). The folder will be created once the main.jl is executed. You may also create it manually.
* ```cesm```: This folder contains the Julia Language source code of the model.
* ```results```: Results will be placed here. The folder will be created when the first results are saved.


# Requirements and Setup
* We use ```Julia 1.6.7```, the long-term-support version. Download it [here](https://julialang.org/downloads/) or [here](https://julialang.org/downloads/oldreleases/).
* We recommend to use a proprietary solver. The model works with [Gurobi](https://www.gurobi.com/solutions/gurobi-optimizer/) and [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer).
* For the initial setup of the CESM, go to the ```cesm``` folder. Select your preferred options in the ```choice_file.jl```, then execute the ```main.jl``` file. All relevant Julia packages then will be instantiated into a seperate Julia environment called ```cesm``` (the folder name) with its files Project.toml and Manifest.toml.
* Only select your preferred features and/or few timesteps. The model otherwise requires prohibitive computation times. 


# Additional data
A simple dataset to test the model is available. Please, ask the authors of the above mentioned publication. Place the files in the ```input_data``` folder in your repo.


# Acknowledgements
The CESM was developped in the project ZellNetz2050 that was funded by the German Federal Ministry for Economic Affairs and Climate Action.  
