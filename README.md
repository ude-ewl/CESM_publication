# Cellular Energy System Model 
This repo contains the Cellular Energy System Model (CESM). The model can be used to perform reserach on market design of decarbonized and decrentralized energy systems.
The CESM model was established for the publication: 
```Schinke-Nendza, A.; Kramer, H.; Khalid, A.; Flatter, F.; Rasti, S.; Uhlemeyer, B; Weber, C. (2024) Modeling of large decentral multi-energy systems: An application to flexible heat pumps in a decarbonized scenario. Journal. DOI:```


# Repo Structure
* ```documentation```: Additional information can be found here 
* ```input_data```: This is the folder of the input data. Data is available at the following data reposity hosted by Zenodo (see below). The folder will be created once the main.jl is executed. You may also create it manually.
* ```model```: This folder contains the Julia Language source code.
* ```results```: Results will be placed here. The folder will be created when the first results are saved.


# Requirements and Setup
* ```Julia 1.6.7``` (long-term-support version). Download the language [here](https://julialang.org/downloads/) or [here](https://julialang.org/downloads/oldreleases/).
* We recommend to use a proprietary solver. The model works with [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer) or [Gurobi](https://www.gurobi.com/solutions/gurobi-optimizer/).
* For the initial setup of the CESM, go to the ```model``` folder. Select options in the ```choice_file.jl```, then execute the ```main.jl``` file. All relevant Julia packages then will be instantiated into a seperate Julia environment called ```model``` with its files Project.toml and Manifest.toml.
* Only select your preferred features and/or few timesteps. The model otherwise requires prohibitive computation times. 


# Additional data
A simple dataset to test the model is available on Zenodo. Please, download the files. Create an ```input_data``` folder in your repo and paste the files there.


# Acknowledgements
The CESM was developped in the project ZellNetz2050 that was funded by the German Federal Ministry for Economic Affairs and Climate Action.  