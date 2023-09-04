module model_chains

using Setfield, SparseArrays

include("market_model.jl")


function calling_solver_package(PARAMETER_SETTINGS)
    if PARAMETER_SETTINGS.SOLVER == "CPLEX"
        eval(:(using CPLEX))
    elseif PARAMETER_SETTINGS.SOLVER == "GLPK"
        eval(:(using GLPK))
    elseif PARAMETER_SETTINGS.SOLVER == "Gurobi"
        eval(:(using Gurobi))
    end
end


function add_constraint_sparse_2D(model, A, x, direction, b, strname)
    
    # obtain index of VariableRef in a very dirty fashion
    string1 = string(b[1,:]) # extract only one dimension (time)
    string2 = replace(string1, r"\n"=>"") # regular expression: remove line breaks
    string3 = replace(string2, r".*:"=>"") # regular expression: remove everything until last occurence of ":" in string
    var_name = replace(split(string3, "[")[1], " "=>"")   
    string4 = replace(replace(string3, (var_name*"[1,")=>""), "]"=>"") # regular expression: remove everything until last occurence of ":" in string
    string5 = split(string4)
    idx = parse.(Int64, string5)
   
    nz_row, nz_col, nz_val = findnz(A) # find nonzero elements, particularly rows to write an constraint for each nonzero row

    collect_vec = collect(1:size(b,1))
    idx_z_row = .!(in.(collect(1:size(b,1)), (nz_row,) )) # if there are empty rows in the sparse assignment matrix
    z_row = collect_vec[idx_z_row]

    if direction == "<="
        @constraint(model, eq_version_nz[unique_nz_row=unique(nz_row),t=idx], sum(nz_val[nz_row.==unique_nz_row] .* x[nz_col[nz_row.==unique_nz_row],t]) <= b[unique_nz_row,t])
    elseif direction == ">="
        @constraint(model, eq_version_nz[unique_nz_row=unique(nz_row),t=idx], sum(nz_val[nz_row.==unique_nz_row] .* x[nz_col[nz_row.==unique_nz_row],t]) >= b[unique_nz_row,t])
    elseif direction == "=="       
        @constraint(model, eq_version_nz[unique_nz_row=unique(nz_row),t=idx], sum(nz_val[nz_row.==unique_nz_row] .* x[nz_col[nz_row.==unique_nz_row],t]) == b[unique_nz_row,t])
    end
    
       if strname != nothing # rename constraints one by one, as general makro doesn't exist in JuMP yet
        counter = 1
        for ct_constraint = unique(nz_row)
            for ct_dim2 = idx
                set_name(eq_version_nz[ct_constraint,ct_dim2], (strname * "[" .*string(counter) * "," * string(ct_dim2) * "]"))
            end
            counter = counter + 1
        end
    end

    # change registration of constraint in model dictionary and unregister generic name to avoid "constraint already defined" error when function is called again
    model[:(strname)] = object_dictionary(model)[:eq_version_nz]
    unregister(model, :eq_version_nz)

    # force nodel variables (b) in rows where there are no units in Ax (and thus there is no enty in sparse matrix) to zero (RHS)
    @constraint(model, eq_version_z_zero[unique_z_row=unique(z_row),t=idx], b[unique_z_row,t] == 0)


    if strname != nothing # rename constraints one by one, as general makro doesn't exist in JuMP yet
        counter = 1
        for ct_constraint = unique(z_row)
            for ct_dim2 = idx
                set_name(eq_version_z_zero[ct_constraint,ct_dim2], (strname * "_zero[" .*string(counter) * "," * string(ct_dim2) * "]"))
            end
            counter = counter + 1
        end
    end

     # change registration of constraint in model dictionary and unregister generic name to avoid "constraint already defined" error when function is called again
     model[:(strname)] = object_dictionary(model)[:eq_version_z_zero]
     unregister(model, :eq_version_z_zero)
 end



#--- Function for rolling optimization
function rolling_optimization(MODEL_INPUT, PARAMETER_SETTINGS)

    # Introduce short notation pointer to same object (as f_model_input is immutable)
    MI = MODEL_INPUT;
    
    # preallocation of "RESULS_ALL" named tuple
    RESULTS_ALL = data_containers.initialize_market_results_all(MI)

    # preallocation of single optimization result variables, as introducing them within the loop would only create local variables that cannot be used in the next iteration
    # should be updated: Save local result to RESULTS_ALL and access data from that variable 
    # pay attention: also required to have a variable to return in the end of the function
    MARKET_RESULT = []

    println("Amount of optimization days: " * string(MI.PERIODS.N_PERIODS))

    #-- Loop for rolling optimization
    for idx_period in 1 : MI.PERIODS.N_PERIODS

        # Check if current period is first period
        IS_FIRST_PERIOD = idx_period == 1 ? 1 : 0

        # Check if current period is last period
        IS_LAST_PERIOD = (idx_period == MI.PERIODS.N_PERIODS) ? 1 : 0

        # Define considered timesteps
        IN_TIME_NO_OVERLAP = ((idx_period - 1) * PARAMETER_SETTINGS.PERIOD_LENGTH + 1) : (idx_period * PARAMETER_SETTINGS.PERIOD_LENGTH)

        # Regularly assign full overlap length, but if close to end, shorten it down
        OVERLAP = PARAMETER_SETTINGS.PERIOD_OVERLAP_LENGTH - (MI.PERIODS.N_PERIODS - idx_period)*PARAMETER_SETTINGS.PERIOD_LENGTH < 0 ? PARAMETER_SETTINGS.PERIOD_OVERLAP_LENGTH : (MI.PERIODS.N_PERIODS - idx_period)*PARAMETER_SETTINGS.PERIOD_LENGTH

        if (IS_LAST_PERIOD == false)
            IN_TIME = ((idx_period - 1) * PARAMETER_SETTINGS.PERIOD_LENGTH + 1) : (idx_period * PARAMETER_SETTINGS.PERIOD_LENGTH + OVERLAP)
        else
            IN_TIME = ((idx_period - 1) * PARAMETER_SETTINGS.PERIOD_LENGTH + 1) : (idx_period * PARAMETER_SETTINGS.PERIOD_LENGTH)
        end
        #println(IN_TIME)

        # set filling levels to zero in first period, else take values of previous optimization period
        TMP_FILLING_LEVEL_0 = (idx_period==1) ? zeros(MI.SETS.STORS[end]) : MARKET_RESULT.FILLING_LEVEL[:,PARAMETER_SETTINGS.PERIOD_LENGTH]
        HP_HS_TMP_FILLING_LEVEL_0 = (idx_period==1) ? zeros(MI.SETS.HPS[end]) : MARKET_RESULT.HP_HS_FILLING_LEVEL[:,PARAMETER_SETTINGS.PERIOD_LENGTH]
        
        # power plant considered to be off before first period, starting at 2nd period, roll over on/off values for unit comittment restrictions from period before
        # (this is done this ways to be able consider UC min online questions, but not implemented intertemporally -> computational burden)
        CONV_ONLINE_LASTPERIOD = (idx_period==1) ? zeros(MI.SETS.GENS[end]) : RESULTS_ALL.CONV_ONLINE_LASTPERIOD[:,idx_period-1]
        CHP_ONLINE_LASTPERIOD = (idx_period==1) ? zeros(MI.SETS.CHPS[end]) : RESULTS_ALL.CHP_ONLINE_LASTPERIOD[:,idx_period-1]

        SETS_LOOP = deepcopy(MI.SETS);
        SETS_LOOP = @set SETS_LOOP.TIME = IN_TIME; # Setfield workaround, as struct SET is principally immuable (for increased performance)
        
        # ========== CESM / MARKET DISPATCH ==========

        global  market_model = build_market_model(MODEL_INPUT, SETS_LOOP, TMP_FILLING_LEVEL_0, IS_FIRST_PERIOD, 
                PARAMETER_SETTINGS, HP_HS_TMP_FILLING_LEVEL_0, CONV_ONLINE_LASTPERIOD, CHP_ONLINE_LASTPERIOD)

        MARKET_RESULT = optimize_market_model(market_model, PARAMETER_SETTINGS, SETS_LOOP)


        #print model status
        println("---\n", now(), " : market model - period ", string(idx_period),"/", string(MI.PERIODS.N_PERIODS), " is ", string(MARKET_RESULT.MODELSTATUS), "\n")

        # print slacks
        if ( sum(MARKET_RESULT.GEN_EL_SLACKS_POS) > 0)
            println("Slack: GEN_EL_SLACKS_POS (sum): ", sum(MARKET_RESULT.GEN_EL_SLACKS_POS), ", ")
        end
        if ( sum(MARKET_RESULT.GEN_EL_SLACKS_NEG) > 0)
            println("Slack: GEN_EL_SLACKS_NEG (sum): ", sum(MARKET_RESULT.GEN_EL_SLACKS_NEG), ", ")
        end
        if ( sum(MARKET_RESULT.GEN_SLACKS_HP) > 0)
            println("Slack: GEN_SLACKS_HP (sum): ", sum(MARKET_RESULT.GEN_SLACKS_HP), ", ")
        end
        if ( sum(MARKET_RESULT.SLACK_Q) > 0)
            println("Slack: SLACK_Q (sum): ", sum(MARKET_RESULT.SLACK_Q), ", ")
        end

        # save market result to list
        RESULTS_ALL = data_containers.append_single_model_result(RESULTS_ALL, MARKET_RESULT, PARAMETER_SETTINGS, idx_period)

    end

    return RESULTS_ALL;
end

end

