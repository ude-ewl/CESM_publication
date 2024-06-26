module choice_file


function input_parameter_settings()

        week_to_start = 9 # testing weeks are 9 (winter, low residual load), 14 (summer, low residual load), 39 (summer, high residual load), 46 (winter, high residual load)
        period_length = 12
        periods = 2

        PARAMETER_SETTINGS = (
            # SLACK_NODE ->Reference node
            SLACK_NODE = 1,
            # PERIOD_LENGTH -> time steps of a single optimization period (w/o overlapping)
            PERIOD_LENGTH = period_length,
            # PERIOD_OVERLAP_LENGTH -> time steps of overlappting time after a single optimization period (to the next period)
            PERIOD_OVERLAP_LENGTH = 12,
            # REF_VOLTAGE ->
            REF_VOLTAGE = 380,
            # CALC_TIME -> time which is considered in the optimization
            CALC_TIME = collect((24*7*(week_to_start-1))+1 : ((24*7*(week_to_start-1))+periods*period_length)),
            # BATTERY STORAGE ON THESE VOTLAGE LEVELS IS MODELLED DETIALED, OTHERS ARE NOT CONSIDERED (PUMPED HYDRO IS ALWAYS CONSIDERED)
            FLEXIBILITY_STORAGE_VOLTAGE_LEVELS = [1], #[1,2,3,4,5], 1=380kV, 2=220kV, 3=110kV, 4=50kV, 5=0.4kV
            # HVDC Losses
            HVDC_LOSSES = "NOT ACTIVE", # "ACTIVE", or "NOT ACTIVE",
            # AC Losses
            AC_LOSSES = "ACTIVE", # "ACTIVE", or "NOT ACTIVE",
            # TRANSFORMER_RESTRICTIONS
            TRANSFORMER_RESTRICTIONS = "NOT ACTIVE", # "ACTIVE", or "NOT ACTIVE",
            # LINE_RESTRICTIONS
            LINE_RESTRICTIONS_VOLTAGE = 110, # exclude line restrictions for lines LOWER than this voltage level, e.g. 0 means that none are excluded (the voltage level written here is the least voltage level included)
            # MARKET DESIGN
            MODEL_TYPE = "CESM", # "CESM" (nodal) or "ZONAL" (uniform price, line limits are ignored)
            # (n-1) thermal restricition of lines (in percent of max thermal line rating)
            N_MINUS_1_FACTOR = 1,
            # (n-1) thermal restricition applies for the line limits of this voltage level and above (the voltage level written here is the least voltage level included)
            N_MINUS_1_VOLTAGE = 220,
            # Choose Solver
            SOLVER = "CPLEX", #"Gurobi", # "CPLEX", "Gurobi" or "GLPK"
            # solver presolve
            PRESOLVE = 1,
            # Flexibility settings for market run
            FLEX_MARKET = (PV_RES="TS", # "Flex"
                            EMOB="TS", 
                            HP="TS"),
                            # STOR is always "Flex" in market run
                            # CHP is always "Flex" in market run
                            # PTG is always "Flex" in market run
                            # PTH is always "Flex" in market run
            
            # Penalty for last timestep to avoid storage emptying
            PENALTY_LAST_TIMESTEP = 5, # EUR/MWh

        )
    return PARAMETER_SETTINGS
end

end