module data_containers

using JuMP

# ===== ONE MARKET RESULTS =====


struct single_model_result
    MODELSTATUS #(1)
    OBJ_VAL
    SOLVE_TIME
    PRICES
    FILLING_LEVEL #(5)
    STORAGES_CHARGING
    STORAGES_DISCHARGING
    OPERATIONAL_COSTS
    GEN_EL_ONLY
    GEN_EL_CHP #(10)
    GEN_HEAT_CHP
    OPERATIONAL_COSTS_CHP
    FLOW_EL_BRANCH
    HVDC_FLOW
    SLACK_Q #(15)
    GEN_EL_BIO
    GEN_HY_ROR
    LOAD_EMOB_RESIDENTIAL
    LOAD_EMOB_COMMERCIAL
    LOAD_EMOB_QUICKCHARGE #(20)
    RENEWABLE_FEED_IN
    GEN_EL_SLACKS_POS
    GEN_EL_SLACKS_NEG
    GEN_PTG
    GEN_HYDROGEN
    GEN_HEAT_PTH
    GEN_HEAT_ONLY
    GEN_HEAT_PUMP
    HP_HS_FILLING_LEVEL
    GENS_ONLINE
    CHPS_ONLINE
    GEN_SLACKS_HP
    HP_P
    HP_HE_P
    HP_P_NODE
    SP_HS_INIT
    SP_HS_WATER_VALUE
    SP_HP_P_MAX
    SP_HE_P_MAX
    SP_HS_Q_MAX_INIT
    SP_HS_Q_MAX
    INVEST_GAS_DUAL #dual of gas investment capacity
    INVEST_PTG_DUAL #dual of ptg investment capacity
end


function single_model_result(PARAMETER_SETTINGS, SETS, model)
    
    if ((PARAMETER_SETTINGS.FLEX_MARKET.HP == "TS") | (PARAMETER_SETTINGS.FLEX_RD.HP == "TS"))
        HP_HS_FILLING_LEVEL = zeros(SETS.HPS[end],PARAMETER_SETTINGS.PERIOD_LENGTH+1)
        HP_P = zeros(SETS.HPS[end],PARAMETER_SETTINGS.PERIOD_LENGTH)
        HP_HE_P = zeros(SETS.HPS[end],PARAMETER_SETTINGS.PERIOD_LENGTH)
        HP_P_NODE = zeros(SETS.NODES[end],PARAMETER_SETTINGS.PERIOD_LENGTH)
        SP_HS_INIT = zeros(SETS.HPS[end],1)
        SP_HS_WATER_VALUE = zeros(SETS.HPS[end],PARAMETER_SETTINGS.PERIOD_LENGTH)
        SP_HP_P_MAX = zeros(SETS.HPS[end],PARAMETER_SETTINGS.PERIOD_LENGTH)
        SP_HE_P_MAX = zeros(SETS.HPS[end],PARAMETER_SETTINGS.PERIOD_LENGTH)
        SP_HS_Q_MAX = zeros(SETS.HPS[end],PARAMETER_SETTINGS.PERIOD_LENGTH+1) # HEAT FILLING LEVEL (SOC)

    elseif ((PARAMETER_SETTINGS.FLEX_MARKET.HP == "Flex"))
        # HP_HS_FILLING_LEVEL = value.(q_hp_hs).data
        HP_HS_FILLING_LEVEL = value.(model[:q_hp_hs]).data
        HP_P = value.(model[:p_hp_unit]).data
        HP_HE_P = value.(model[:p_hp_he_unit]).data
        HP_P_NODE = value.(model[:p_hp_node]).data
        SP_HS_INIT = dual.(model[:cons_heatpump_heatstorage_initial_filling_level]).data # 1 dual value
        SP_HS_WATER_VALUE = dual.(model[:cons_heatstorage_filling_level]).data           # 24 dual values
        SP_HP_P_MAX = dual.(model[:cons_heatpump_max_power]).data
        SP_HE_P_MAX = dual.(model[:cons_heatpump_heatelement_max_power]).data
        SP_HS_Q_MAX = dual.(model[:cons_heatpump_heatstorage_max_capacity]).data         # 25 dual values
    end

    if PARAMETER_SETTINGS.MODEL_TYPE != "BENDERS INVEST"
        INVEST_GAS_DUAL = zeros(SETS.GENS[end],1)
        INVEST_PTG_DUAL = zeros(SETS.PTG[end],1)
    else
        #println(value.(model[:p_g_unit_invest]).data)
        INVEST_GAS_DUAL = dual.(model[:cons_conventional_invest]).data
        INVEST_PTG_DUAL = dual.(model[:cons_ptg_invest]).data

    end
    #println(PARAMETER_SETTINGS.MODEL_TYPE != "BENDERS INVEST")

    return single_model_result(   
        termination_status(model), #(1)
        objective_value(model),
        solve_time(model),
        dual.(model[:cons_load_serving_elec_export]).data, #prices
        value.(model[:e_s]).data, #(5)
        value.(model[:p_ch]).data, 
        value.(model[:p_dis]).data,
        value.(model[:c_g_op]).data,
        value.(model[:p_g_unit]).data,
        value.(model[:p_chp_unit]).data, #(10)
        value.(model[:q_chp_unit]).data,
        value.(model[:c_chp_op]).data,
        value.(model[:p_branch]).data, # branch hvac flow
        value.(model[:p_hvdc_line]).data, # branch hvdc flow
        value.(model[:q_chp_region_slack]).data, #(15)
        value.(model[:p_bio_unit]).data,
        value.(model[:p_hy_ror_unit]).data,
        value.(model[:p_emob_node_residential]).data,
        value.(model[:p_emob_node_commercial]).data, 
        value.(model[:p_emob_node_quick]).data, #(20)
        value.(model[:p_re]).data,
        value.(model[:p_slack_pos]).data, #(25)
        value.(model[:p_slack_neg]).data,
        value.(model[:p_ptg_unit]).data,
        value.(model[:p_ptg_unit_hydrogen]).data,
        value.(model[:p_pth_unit]).data, #(30)
        value.(model[:p_heat_only_unit]).data,
        value.(model[:p_hp_node]).data,
		HP_HS_FILLING_LEVEL,
        value.(model[:p_g_unit_online]).data,
        value.(model[:p_chp_unit_online]).data,
        value.(model[:hp_heat_slack]).data,
        HP_P, #(35)
        HP_HE_P,
        HP_P_NODE,
        SP_HS_INIT,
        SP_HS_WATER_VALUE,
        SP_HP_P_MAX, #(40)
        SP_HE_P_MAX,
        SP_HS_Q_MAX[:,1], # init value
        SP_HS_Q_MAX[:,2:end],
        INVEST_GAS_DUAL, # dual of invested gas capacity (relevant for Benders decomposition)
        INVEST_PTG_DUAL # dual of invested ptg capacity (relevant for Benders decomposition)
        )
end



# ===== ALL MARKET RESULTS -- INITIALIZATION =====

function initialize_market_results_all(MI)

    market_keys = ( :MODELSTATUS, #(1)
        :OBJ_VAL, 
        :SOLVE_TIME,
        :PRICES,
        :CONV_ONLINE_LASTPERIOD, #(5)
        :CHP_ONLINE_LASTPERIOD,
        :TMP_FILLING_LEVEL_0,
        :OPERATIONAL_COSTS,
        :FILLING_LEVEL,
        :GEN_EL_ONLY, #(10)
        :GEN_EL_CHP,
        :GEN_HEAT_CHP,
        :OPERATIONAL_COSTS_CHP,
        :FLOW_EL_BRANCH,
        :HVDC_FLOW, #(15)
        :SLACK_Q,
        :GEN_EL_BIO,
        :GEN_HY_ROR,
        :LOAD_EMOB_RESIDENTIAL,
        :LOAD_EMOB_COMMERCIAL, #(20)
        :LOAD_EMOB_QUICKCHARGE,
        :RENEWABLE_FEED_IN,
        :STORAGES_CHARGING,
        :STORAGES_DISCHARGING,
        :GEN_EL_SLACKS_POS, #(25)
        :GEN_EL_SLACKS_NEG,
        :HP_HS_TMP_FILLING_LEVEL_0,
        :GEN_PTG,
        :GEN_HYDROGEN,
        :GEN_HEAT_PTH, #(30)
        :GEN_HEAT_ONLY,
        :GEN_HEAT_PUMP,
        :GEN_SLACKS_HP,
        :HP_P,
        :HP_HE_P, #(35)
        :HP_P_NODE,
        :HP_HS_FILLING_LEVEL,
        :SP_HS_INIT,
        :SP_HS_WATER_VALUE,
        :SP_HP_P_MAX, #(40)
        :SP_HE_P_MAX,
        :SP_HS_Q_MAX_INIT,
        :SP_HS_Q_MAX,
        :INVEST_GAS_DUAL, # dual of invested gas capacity (relevant for Benders decomposition)
        :INVEST_PTG_DUAL # dual of invested ptg capacity (relevant for Benders decomposition)
    )

    empty_vals = ( Vector{String}(broadcast(string, zeros(MI.PERIODS.N_PERIODS))), # MODELSTATUS
        Array{Union{Float64,Missing}}(missing,MI.PERIODS.N_PERIODS), # OBJ_VAL
        Array{Union{Float64,Missing}}(missing,MI.PERIODS.N_PERIODS), # SOLVE TIME
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # PRICES
        Array{Union{Float64,Missing}}(missing,MI.SETS.GENS[end], MI.PERIODS.N_PERIODS), # CONV_ONLINE_LASTPERIOD
        Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.PERIODS.N_PERIODS), # CHP_ONLINE_LASTPERIOD
        Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]+1), # TMP_FILLING_LEVEL_0
        Array{Union{Float64,Missing}}(missing,MI.SETS.GENS[end], MI.SETS.TIME[end]),  # OPERATIONAL_COSTS
        Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]), # FILLING_LEVEL
        Array{Union{Float64,Missing}}(missing,MI.SETS.GENS[end], MI.SETS.TIME[end]), # GEN_EL_ONLY
        Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.SETS.TIME[end]), # GEN_EL_CHP
        Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.SETS.TIME[end]), # GEN_HEAT_CHP
        Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.SETS.TIME[end]), # CHP Operational Costs
        Array{Union{Float64,Missing}}(missing,MI.SETS.BRANCHES[end], MI.SETS.TIME[end]), # FLOW_EL_BRANCH
        Array{Union{Float64,Missing}}(missing,MI.SETS.HVDCS[end], MI.SETS.TIME[end]), # HVDC_FLOW
        Array{Union{Float64,Missing}}(missing,MI.SETS.HEATREGIONS[end], MI.SETS.TIME[end]), #  SLACK HEAT
        Array{Union{Float64,Missing}}(missing,MI.SETS.BIOS[end], MI.SETS.TIME[end]), # GEN_EL_BIO
        Array{Union{Float64,Missing}}(missing,MI.SETS.HY_ROR[end], MI.SETS.TIME[end]), # GEN_HY_ROR
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # EMOB RESIDENTIAL at NODE
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # EMOB COMMERCIAL at NODE
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # EMOB QUICK at NODE
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # RENEWABLE_FEED_IN
        Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]), # STORAGES_CHARGING
        Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]), # STORAGES_DISCHARGING
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # GEN_EL_SLACKS_POS
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # GEN_EL_SLACKS_NEG
        # zeros(SETS.HPS[end]), # HP_HS_TMP_FILLING_LEVEL_0
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]+1), # HP_HS_TMP_FILLING_LEVEL_0
        Array{Union{Float64,Missing}}(missing,MI.SETS.PTG[end], MI.SETS.TIME[end]), # GEN_PTG
        Array{Union{Float64,Missing}}(missing,MI.SETS.PTG[end], MI.SETS.TIME[end]), # GEN_HYDROGEN
        Array{Union{Float64,Missing}}(missing,MI.SETS.PTH[end], MI.SETS.TIME[end]), # GEN_HEAT_PTH
        Array{Union{Float64,Missing}}(missing,MI.SETS.HEAT_ONLY[end], MI.SETS.TIME[end]), # GEN_HEAT_ONLY
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), #GEN_HEAT_PUMP
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), # GEN_SLACKS_HP
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #HP_P
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #HP_HE_P
        Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # HP_P_NODE
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #HP_HS_FILLING_LEVEL
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.PERIODS.N_PERIODS), # SP_HS_INIT
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), # SP_HS_WATER_VALUE
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #SP_HP_P_MAX
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), # SP_HE_P_MAX
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.PERIODS.N_PERIODS), # SP_HS_Q_MAX_INIT
        Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), # SP_HS_Q_MAX
        Array{Union{Float64,Missing}}(missing,MI.SETS.GENS[end]), # INVEST_GAS_DUAL
        Array{Union{Float64,Missing}}(missing,MI.SETS.PTH[end]) # INVEST_PTG_DUAL
    )

    return (; zip(market_keys, empty_vals)...)
end

function initialize_rd_results_all(MI)

    rd_keys = ( 
        :MODELSTATUS, #(1)
        :OBJ_VAL, 
        :SOLVE_TIME,
        :PRICES,
        :CONV_ONLINE_LASTPERIOD, #(5)
        :CHP_ONLINE_LASTPERIOD,
        :TMP_FILLING_LEVEL_0,
        :OPERATIONAL_COSTS,
        :FILLING_LEVEL,
        :GEN_EL_ONLY, #(10)
        :GEN_EL_CHP,
        :GEN_HEAT_CHP,
        :OPERATIONAL_COSTS_CHP,
        :FLOW_EL_BRANCH,
        :HVDC_FLOW, #(15)
        :SLACK_Q,
        :GEN_EL_BIO,
        :GEN_HY_ROR,
        :LOAD_EMOB_RESIDENTIAL,
        :LOAD_EMOB_COMMERCIAL, #(20)
        :LOAD_EMOB_QUICKCHARGE,
        :RENEWABLE_FEED_IN,
        :STORAGES_CHARGING,
        :STORAGES_DISCHARGING,
        :GEN_EL_SLACKS_POS, #(25)
        :GEN_EL_SLACKS_NEG,
        :HP_HS_TMP_FILLING_LEVEL_0,
        :GEN_PTG,
        :GEN_HYDROGEN,
        :GEN_HEAT_PTH, #(30)
        :GEN_HEAT_ONLY,
        :GEN_HEAT_PUMP,
        :GEN_SLACKS_HP,
        :HP_P,
        :HP_HE_P, #(35)
        :HP_P_NODE,
        :HP_HS_FILLING_LEVEL,
        :SP_HS_INIT,
        :SP_HS_WATER_VALUE,
        :SP_HP_P_MAX, #(40)
        :SP_HE_P_MAX,
        :SP_HS_Q_MAX_INIT,
        :SP_HS_Q_MAX
        )

    empty_vals = ( Vector{String}(broadcast(string, zeros(MI.PERIODS.N_PERIODS))), # MODELSTATUS
    Array{Union{Float64,Missing}}(missing,MI.PERIODS.N_PERIODS), # OBJ_VAL
    Array{Union{Float64,Missing}}(missing,MI.PERIODS.N_PERIODS), # SOLVE TIME
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # PRICES
    Array{Union{Float64,Missing}}(missing,MI.SETS.GENS[end], MI.PERIODS.N_PERIODS), # CONV_ONLINE_LASTPERIOD
    Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.PERIODS.N_PERIODS), # CHP_ONLINE_LASTPERIOD
    Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]+1), # TMP_FILLING_LEVEL_0
    Array{Union{Float64,Missing}}(missing,MI.SETS.GENS[end], MI.SETS.TIME[end]),  # OPERATIONAL_COSTS
    Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]), # FILLING_LEVEL
    Array{Union{Float64,Missing}}(missing,MI.SETS.GENS[end], MI.SETS.TIME[end]), # GEN_EL_ONLY
    Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.SETS.TIME[end]), # GEN_EL_CHP
    Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.SETS.TIME[end]), # GEN_HEAT_CHP
    Array{Union{Float64,Missing}}(missing,MI.SETS.CHPS[end], MI.SETS.TIME[end]), # CHP Operational Costs
    Array{Union{Float64,Missing}}(missing,MI.SETS.BRANCHES[end], MI.SETS.TIME[end]), # FLOW_EL_BRANCH
    Array{Union{Float64,Missing}}(missing,MI.SETS.HVDCS[end], MI.SETS.TIME[end]), # HVDC_FLOW
    Array{Union{Float64,Missing}}(missing,MI.SETS.HEATREGIONS[end], MI.SETS.TIME[end]), #  SLACK HEAT
    Array{Union{Float64,Missing}}(missing,MI.SETS.BIOS[end], MI.SETS.TIME[end]), # GEN_EL_BIO
    Array{Union{Float64,Missing}}(missing,MI.SETS.HY_ROR[end], MI.SETS.TIME[end]), # GEN_HY_ROR
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # EMOB RESIDENTIAL at NODE
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # EMOB COMMERCIAL at NODE
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # EMOB QUICK at NODE
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # RENEWABLE_FEED_IN
    Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]), # STORAGES_CHARGING
    Array{Union{Float64,Missing}}(missing,MI.SETS.STORS[end], MI.SETS.TIME[end]), # STORAGES_DISCHARGING
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # GEN_EL_SLACKS_POS
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # GEN_EL_SLACKS_NEG
    # zeros(SETS.HPS[end]), # HP_HS_TMP_FILLING_LEVEL_0
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]+1), # HP_HS_TMP_FILLING_LEVEL_0
    Array{Union{Float64,Missing}}(missing,MI.SETS.PTG[end], MI.SETS.TIME[end]), # GEN_PTG
    Array{Union{Float64,Missing}}(missing,MI.SETS.PTG[end], MI.SETS.TIME[end]), # GEN_HYDROGEN
    Array{Union{Float64,Missing}}(missing,MI.SETS.PTH[end], MI.SETS.TIME[end]), # GEN_HEAT_PTH
    Array{Union{Float64,Missing}}(missing,MI.SETS.HEAT_ONLY[end], MI.SETS.TIME[end]), # GEN_HEAT_ONLY
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), #GEN_HEAT_PUMP
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), # GEN_SLACKS_HP
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #HP_P
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #HP_HE_P
    Array{Union{Float64,Missing}}(missing,MI.SETS.NODES[end], MI.SETS.TIME[end]), # HP_P_NODE
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #HP_HS_FILLING_LEVEL
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.PERIODS.N_PERIODS), # SP_HS_INIT
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), # SP_HS_WATER_VALUE
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), #SP_HP_P_MAX
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]), # SP_HE_P_MAX
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.PERIODS.N_PERIODS), # SP_HS_Q_MAX_INIT
    Array{Union{Float64,Missing}}(missing,MI.SETS.HPS[end], MI.SETS.TIME[end]) # SP_HS_Q_MAX
    )

    return (; zip(rd_keys, empty_vals)...)
end


# ===== ALL MARKET RESULTS -- FILL WITH ITERATION DATA =====
function append_single_model_result(RESULTS_ALL, MARKET_RESULT, PARAMETER_SETTINGS, idx_period)

    IN_TIME_NO_OVERLAP = ((idx_period - 1) * PARAMETER_SETTINGS.PERIOD_LENGTH + 1) : (idx_period * PARAMETER_SETTINGS.PERIOD_LENGTH)

         # Save results of optimization: idx_period
         RESULTS_ALL.MODELSTATUS[idx_period] = string(MARKET_RESULT.MODELSTATUS)
         RESULTS_ALL.OBJ_VAL[idx_period] = MARKET_RESULT.OBJ_VAL
         RESULTS_ALL.SOLVE_TIME[idx_period] = MARKET_RESULT.SOLVE_TIME
         RESULTS_ALL.PRICES[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.PRICES[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.CONV_ONLINE_LASTPERIOD[:, idx_period] = MARKET_RESULT.GENS_ONLINE[:,PARAMETER_SETTINGS.PERIOD_LENGTH] # only last value
         RESULTS_ALL.CHP_ONLINE_LASTPERIOD[:, idx_period] = MARKET_RESULT.CHPS_ONLINE[:,PARAMETER_SETTINGS.PERIOD_LENGTH] # only last value
         RESULTS_ALL.SP_HS_INIT[:, idx_period] = MARKET_RESULT.SP_HS_INIT
         RESULTS_ALL.OPERATIONAL_COSTS[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.OPERATIONAL_COSTS[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.FILLING_LEVEL[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.FILLING_LEVEL[:,1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         
         # Save results of optimization: IN_TIME_NO_OVERLAP
         RESULTS_ALL.GEN_EL_ONLY[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_EL_ONLY[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_EL_CHP[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_EL_CHP[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_HEAT_CHP[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_HEAT_CHP[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.OPERATIONAL_COSTS_CHP[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.OPERATIONAL_COSTS_CHP[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.FLOW_EL_BRANCH[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.FLOW_EL_BRANCH[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.HVDC_FLOW[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.HVDC_FLOW[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.SLACK_Q[:,IN_TIME_NO_OVERLAP] = MARKET_RESULT.SLACK_Q[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_EL_BIO[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_EL_BIO[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_HY_ROR[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_HY_ROR[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.LOAD_EMOB_RESIDENTIAL[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.LOAD_EMOB_RESIDENTIAL[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.LOAD_EMOB_COMMERCIAL[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.LOAD_EMOB_COMMERCIAL[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.LOAD_EMOB_QUICKCHARGE[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.LOAD_EMOB_QUICKCHARGE[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.RENEWABLE_FEED_IN[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.RENEWABLE_FEED_IN[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.STORAGES_CHARGING[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.STORAGES_CHARGING[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.STORAGES_DISCHARGING[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.STORAGES_DISCHARGING[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_EL_SLACKS_POS[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_EL_SLACKS_POS[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_EL_SLACKS_NEG[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_EL_SLACKS_NEG[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_PTG[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_PTG[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_HYDROGEN[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_HYDROGEN[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_HEAT_PTH[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_HEAT_PTH[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_HEAT_ONLY[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_HEAT_ONLY[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_HEAT_PUMP[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_HEAT_PUMP[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.GEN_SLACKS_HP[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.GEN_SLACKS_HP[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.HP_P[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.HP_P[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.HP_HE_P[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.HP_HE_P[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.HP_P_NODE[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.HP_P_NODE[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
            RESULTS_ALL.HP_HS_FILLING_LEVEL[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.HP_HS_FILLING_LEVEL[:,1:PARAMETER_SETTINGS.PERIOD_LENGTH]
            RESULTS_ALL.SP_HS_INIT[:, idx_period] = MARKET_RESULT.HP_HS_FILLING_LEVEL[:,1]
         RESULTS_ALL.SP_HS_WATER_VALUE[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.SP_HS_WATER_VALUE[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.SP_HP_P_MAX[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.SP_HP_P_MAX[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.SP_HE_P_MAX[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.SP_HE_P_MAX[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
         RESULTS_ALL.SP_HS_Q_MAX[:, IN_TIME_NO_OVERLAP] = MARKET_RESULT.SP_HS_Q_MAX[:, 1:PARAMETER_SETTINGS.PERIOD_LENGTH]
            RESULTS_ALL.SP_HS_Q_MAX_INIT[:, idx_period] = MARKET_RESULT.SP_HS_Q_MAX_INIT

    return RESULTS_ALL

end





end