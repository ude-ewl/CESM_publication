

using JuMP, XLSX, CSV, DataFrames, LinearAlgebra, Dates, BenchmarkTools, MarketTechnicals, Dates

include("data_containers.jl")

#--- Create optimization function
function build_market_model(MI, SETS, E_0, IS_FIRST_PERIOD, PARAMETER_SETTINGS, HP_Q_0, CONV_ONLINE_LASTPERIOD, CHP_ONLINE_LASTPERIOD)

    c_invest_gas = 999*1e3/52
    c_invest_ptg = 9*1e3/52

    BRANCHES = MI.BRANCHES
    LOADS = MI.LOADS
    RENEWABLES = MI.RENEWABLES
    CONVENTIONALS = MI.CONVENTIONALS
    STORAGES = MI.STORAGES
    ASSIGNMENTS = MI.ASSIGNMENTS
    PTG = MI.PTG
    COMBINEDHEATPOWERS = MI.COMBINEDHEATPOWERS
    HEAT_DEMAND =  MI.HEAT_DEMAND
    HVDCS = MI. HVDCS
    BIOS = MI.BIOS 
    HY_ROR = MI.HY_ROR
    EMOB_RESIDENTIAL = MI.EMOB_RESIDENTIAL  
    EMOB_COMMERCIAL = MI.EMOB_COMMERCIAL
    EMOB_QUICK = MI.EMOB_QUICK
    HEATPUMPS = MI.HEATPUMPS
    PTH = MI.PTH
    HEAT_ONLY = MI.HEAT_ONLY
    CONST_HEAT_LOAD = MI.CONST_HEAT_LOAD

    #--- Initialize optimization
    # Initialize Optimization model
    if PARAMETER_SETTINGS.SOLVER == "CPLEX"
        model = Model(optimizer_with_attributes(CPLEX.Optimizer,
            "CPXPARAM_ScreenOutput" => 0, 
            "CPXPARAM_Simplex_Tolerances_Markowitz" => 0.1,
            "CPXPARAM_Simplex_Tolerances_Optimality" => 1e-2,
            "CPXPARAM_Output_WriteLevel" => 3))
    elseif PARAMETER_SETTINGS.SOLVER == "GLPK"
        model = Model(GLPK.Optimizer)
    elseif PARAMETER_SETTINGS.SOLVER == "Gurobi"
        model = Model(optimizer_with_attributes(Gurobi.Optimizer,
            "output_flag" => false, ))
        set_optimizer_attribute(model, "Presolve", PARAMETER_SETTINGS.PRESOLVE)
        set_optimizer_attribute(model, "Method", 1) # dual simplex
    end


    println("---\n", now(), " : market model: build...")
    #--- Optimization: Variable definitions

    #-- Demand
    # Electricity demand during SETS.TIME step t at node n (in MW)
    @variable(model, p_d[SETS.NODES,SETS.TIME])

    # Electricity export (negative values imply an import) at node n during SETS.TIME step t (in MW)
    @variable(model, p_ex[SETS.NODES,SETS.TIME])
    # Slack for electricity export
    @variable(model, p_slack_pos[SETS.NODES,SETS.TIME]>= 0)
    @variable(model, p_slack_neg[SETS.NODES,SETS.TIME]>= 0)

    #-- Power flow
    # Power flow across BRANCH m ∈ M during SETS.TIME step t (in MW)
    @variable(model, p_branch[SETS.BRANCHES,SETS.TIME])
    # Absolute value of power flow across BRANCH m ∈ M during SETS.TIME step t (in MW)
    @variable(model, 0 <= p_branch_abs[SETS.BRANCHES,SETS.TIME] )
    # voltage angle delta_lf
    @variable(model, delta_lf[SETS.NODES,SETS.TIME])
    if PARAMETER_SETTINGS.AC_LOSSES == "ACTIVE"
        # Electricity losses at node n during TIME step t (in MW)
        @variable(model, 0 <= p_loss_nodal[SETS.NODES,SETS.TIME] )
        @variable(model, 0 <= p_loss_line_abs[m=SETS.BRANCHES,t=SETS.TIME] <= (BRANCHES.LOSS_FACTOR_m[m] * BRANCHES.PMAX[m]) * 2 )
    end

    #-- HVDC
    # HVDC power exchange at designated node
    @variable(model, p_hvdc_node[SETS.NODES,SETS.TIME])
    # HVDC line loading
    @variable(model, p_hvdc_line[SETS.HVDCS,SETS.TIME])
    if PARAMETER_SETTINGS.HVDC_LOSSES == "ACTIVE"
        # Absolute power flow across HVDC line
        @variable(model, p_hvdc_abs[SETS.HVDCS,SETS.TIME])
        # HVDC losses is solely affine dependent on p_hvdc_line
        @variable(model, p_hvdc_loss[SETS.HVDCS,SETS.TIME] >= 0)
    end

    #-- Power to gas (hydrogen)
    # Production of PTG unit k during SETS.TIME step t at node n (in MW)
    @variable(model, p_ptg_unit[k=SETS.PTG, t=SETS.TIME] >= 0)
    # Nodal electricity consumption of PTG units
    @variable(model, p_ptg_node[k=SETS.NODES, t=SETS.TIME] >= 0)
    # Hydrogen production of PTG units
    @variable(model, p_ptg_unit_hydrogen[k=SETS.PTG, t=SETS.TIME] >= 0)
    # Slack generator for power to gas
    @variable(model, p_ptg_slack[t=SETS.TIME] >= 0)
    # # Online capacity of PTG units
    # @variable(model, p_ptg_unit_hydrogen_online[k=SETS.PTG, t=SETS.TIME] >= 0)

    #-- Conventional (electricity only) power plants (generators)
    # Production of power plant k during SETS.TIME step t at node n (in MW)
    @variable(model, p_g_unit[SETS.GENS,SETS.TIME] >= 0)
    # Start up capacity
    @variable(model, p_g_unit_startup[SETS.GENS, SETS.TIME] >= 0)
    # (SETS.TIME[begin]+1):SETS.TIME[end]] >= 0)
    # (SETS.TIME[begin]-1):SETS.TIME[end]] >= 0)
    # [Auxiliary variable] Available (online) power of power plant k during SETS.TIME step t at node n (in MW)
    @variable(model, p_g_unit_online[SETS.GENS,SETS.TIME] >= 0)
    # Nodal power injection of power plants at node n during SETS.TIME step t (in MW)
    @variable(model, p_g_node[SETS.NODES,SETS.TIME] >= 0)
    #- Costs
    # Operational costs of power plant k at SETS.TIME step t (in €)
    @variable(model, c_g_op[SETS.GENS, SETS.TIME] >= 0)

    # investment variable if endogeneous gas capacity modeling
    if ((PARAMETER_SETTINGS.MODEL_TYPE == "INVEST"))
        @variable(model, p_g_unit_invest[SETS.GENS] >= 0)
        @variable(model, p_ptg_unit_invest[k=SETS.PTG] >= 0)
    end


    #-- Renewables
    # Nodal power injection of renewables at node n during SETS.TIME step t (in MW)
    @variable(model, p_re[SETS.NODES,SETS.TIME] >= 0)
    @variable(model, p_pv_wind_nodal[SETS.NODES,SETS.TIME] >= 0)
    @variable(model, p_pv_wind_after_curtailment[SETS.RENS,SETS.TIME] >= 0)

    #-- Hydro run-of-river
    # Nodal power injection of power plants at node n during SETS.TIME step t (in MW)
    @variable(model, p_hy_ror_nodal[SETS.NODES, SETS.TIME] >= 0)
    # # Production of power plant k during SETS.TIME step t at node n (in MW)
    @variable(model, p_hy_ror_unit[SETS.HY_ROR, SETS.TIME] >= 0)

    #-- Combined heat power plants (CHP)
    # Electricity production of CHP k during SETS.TIME step t at node n (in MW)
    @variable(model, p_chp_unit[SETS.CHPS,SETS.TIME] >= 0)
    # Start up capacity
    @variable(model, p_chp_unit_startup[SETS.CHPS,SETS.TIME] >= 0)
    # (SETS.TIME[begin]+1):SETS.TIME[end]] >= 0)
    # (SETS.TIME[begin]-1):SETS.TIME[end]] >= 0)
    # Nodal power injection of CHPs at node n during SETS.TIME step t (in MW)
    @variable(model, p_chp_node[SETS.NODES,SETS.TIME] >= 0)
    # Heat production of CHP k during SETS.TIME step t at node n (in MW)
    @variable(model, q_chp_unit[SETS.CHPS,SETS.TIME] >= 0)
    # Nodal heat injection of CHPs in heat region r during SETS.TIME step t (in MW)
    @variable(model, q_chp_region[SETS.HEATREGIONS,SETS.TIME] >= 0)
    # Nodal slack heat injection of CHPs in heat region r during SETS.TIME step t (in MW)
    @variable(model, q_chp_region_slack[SETS.HEATREGIONS,SETS.TIME] >= 0)
    # [Auxiliary variable] Available (online) power of CHP k during SETS.TIME step t at node n (in MW)
    @variable(model, p_chp_unit_online[SETS.CHPS,SETS.TIME] >= 0)
    #- Costs
    # Operational costs of CHP k at SETS.TIME step t (in €)
    @variable(model, c_chp_op[SETS.CHPS,SETS.TIME] >= 0)

    #-- Heat demand
    # Nodal heat demand during SETS.TIME step t at node n (in MW)
    @variable(model, q_heat_demand[SETS.HEATREGIONS,SETS.TIME] >= 0)

    # -- Power to heat
    # Heat prodution of units
    @variable(model, p_pth_unit[SETS.PTH, SETS.TIME] >= 0)
    # Heat prodution in heat regions
    @variable(model, p_pth_region[SETS.HEATREGIONS, SETS.TIME] >= 0)
    # Electricty consumption of units for elec_nodes
    @variable(model, p_pth_node[SETS.NODES, SETS.TIME] >= 0)

    # -- Heat only
    # Heat production of units
    @variable(model, p_heat_only_unit[SETS.HEAT_ONLY, SETS.TIME] >= 0)
    # Heat prodution in heat regions
    @variable(model, p_heat_only_region[SETS.HEATREGIONS, SETS.TIME] >= 0)

    # -- Constant heat load
    # Constant heat load/demand for heat region
    @variable(model, p_const_heat_load_region[SETS.HEATREGIONS, SETS.TIME] >= 0)

    #-- Biomass/gas
    # Production of power plant k during SETS.TIME step t at node n (in MW)
    @variable(model, p_bio_unit[SETS.BIOS, SETS.TIME] >= 0)
    # Nodal power injection of power plants at node n during SETS.TIME step t (in MW)
    @variable(model, p_bio_node[SETS.NODES, SETS.TIME] >= 0)
    # Costs
    @variable(model, c_bio_op[SETS.BIOS, SETS.TIME] >= 0)





    #-- Heat pumps
    # Comsumption of heat pumps at node n during timestep t (in MW)
    @variable(model, p_hp_node[SETS.NODES, SETS.TIME] >= 0)
    # If detailed model is selected for flexibilities
    if PARAMETER_SETTINGS.FLEX_MARKET.HP == "Flex"
        # Electricuty consumption of heat pump s during timestep t (in MW)
        @variable(model, p_hp_unit[s=SETS.HPS, t=SETS.TIME] >= 0)
        # Electricty consumption of heat element s during timestep t (in MW)
        @variable(model, p_hp_he_unit[s=SETS.HPS, t=SETS.TIME] >= 0)
        # both together
        @variable(model, p_hp_and_he_unit[s=SETS.HPS, t=SETS.TIME] >= 0)
        # COP of heat pump s during timestep t (in MW)
        @variable(model, cop_hp[s=SETS.HPS, t=SETS.TIME])
        # Heat demand of heat pump s during timestep t (in MW)
        @variable(model, q_hp_d[s=SETS.HPS, t=SETS.TIME])
        # Filling level of storage k during SETS.TIME step t at node n (in MW * t)
        @variable(model, q_hp_hs[s=SETS.HPS, (SETS.TIME[begin]-1):SETS.TIME[end]] >= 0)
    end
    # Heat pump slack
    @variable(model, hp_heat_slack[s=SETS.HPS, t=SETS.TIME] >= 0)

    #-- Storage units
    # Nodal power injection of storages at node n during SETS.TIME step t (in MW)
    @variable(model, p_s[SETS.NODES,SETS.TIME])
    # Charge of storage k during SETS.TIME step t at unit k (in MW)
    @variable(model, 0 <= p_ch[SETS.STORS,SETS.TIME] )
    # Discharge of storage k during SETS.TIME step t at unit k (in MW)
    @variable(model, 0 <= p_dis[SETS.STORS,SETS.TIME])
    # Storage power exchange at unit k (in MW)
    @variable(model, p_dis_minus_ch[k=SETS.STORS,t=SETS.TIME])
    # Filling level of storage k during SETS.TIME step t at node n (in MW * t)
    @variable(model, e_s[SETS.STORS,(SETS.TIME[begin]-1):SETS.TIME[end]] >= 0)

    #-- Electric mobility
    # Comsuption of electic mobility at node n during SETS.TIME step t (in MW)
    @variable(model, p_emob_node[SETS.NODES, SETS.TIME] >= 0)
    @variable(model, p_emob_node_residential[SETS.NODES, SETS.TIME] >= 0)
    @variable(model, p_emob_node_commercial[SETS.NODES, SETS.TIME] >= 0)
    @variable(model, p_emob_node_quick[SETS.NODES, SETS.TIME] >= 0)

    if PARAMETER_SETTINGS.FLEX_MARKET.EMOB == "Flex"
        @variable(model,p_emob_node_residential_max[SETS.NODES, SETS.TIME] >= 0)
        @variable(model,p_emob_node_residential_dumb_dem[SETS.NODES, SETS.TIME] >= 0)
    end

    #--- Optimization: Constraint definitions

    if PARAMETER_SETTINGS.MODEL_TYPE == "CESM" # consider line restrictions if CESM (nodal) setup, when zonal+RD, do NOT use line restrictions
        # Get absolute power flow across branches
        @constraint(model, cons_line_losses_abs_1[m=SETS.BRANCHES, t=SETS.TIME], p_branch_abs[m,t] >= p_branch[m,t] )
        @constraint(model, cons_line_losses_abs_2[m=SETS.BRANCHES, t=SETS.TIME], p_branch_abs[m,t] >= -p_branch[m,t] )
    end
    #-- Losses
    if PARAMETER_SETTINGS.AC_LOSSES == "ACTIVE" && PARAMETER_SETTINGS.MODEL_TYPE == "CESM"
        # Do not consider losses in case of a zonal run
        @constraint(model, cons_line_losses[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] == BRANCHES.LOSS_FACTOR_m[m] * p_branch_abs[m,t] )
        # # # @constraint(model, cons_line_losses_1[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR_m0[m] * p_branch_abs[m,t] + BRANCHES.LOSS_FACTOR_b0[m] )
        # # @constraint(model, cons_line_losses_1[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR_m1[m] * p_branch_abs[m,t] + BRANCHES.LOSS_FACTOR_b1[m] )
        # # @constraint(model, cons_line_losses_2[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR_m2[m] * p_branch_abs[m,t] + BRANCHES.LOSS_FACTOR_b2[m] )
        # # @constraint(model, cons_line_losses_3[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR_m3[m] * p_branch_abs[m,t] + BRANCHES.LOSS_FACTOR_b3[m] )
        # # @constraint(model, cons_line_losses_4[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR_m4[m] * p_branch_abs[m,t] + BRANCHES.LOSS_FACTOR_b4[m] )
        # @constraint(model, cons_line_losses_1[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR_m1a[m] * p_branch_abs[m,t] + BRANCHES.LOSS_FACTOR_b1a[m] )
        # @constraint(model, cons_line_losses_2[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR_m2a[m] * p_branch_abs[m,t] + BRANCHES.LOSS_FACTOR_b2a[m] )
        # # # @constraint(model, cons_line_losses_abs_1[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= BRANCHES.LOSS_FACTOR[m] .* p_branch[m,t] )
        # # # @constraint(model, cons_line_losses_abs_2[m=SETS.BRANCHES, t=SETS.TIME], p_loss_line_abs[m,t] >= -BRANCHES.LOSS_FACTOR[m] .* p_branch[m,t] )
        @constraint(model, cons_line_losses_nodal[n=SETS.NODES, t=SETS.TIME], p_loss_nodal[n,t] == ((BRANCHES.ASSIGNMENT_FROM[n,:] + 
            BRANCHES.ASSIGNMENT_TO[n,:])./2)' * p_loss_line_abs[:,t] )
        @expression(model, expression_line_losses_nodal[n=SETS.NODES, t=SETS.TIME], p_loss_nodal[n,t] )
        # @expression(model, expression_line_losses_nodal[n=SETS.NODES, t=SETS.TIME], ((BRANCHES.ASSIGNMENT_FROM[n,:] + 
        #     BRANCHES.ASSIGNMENT_TO[n,:])./2)' * p_loss_line_abs[:,t] )
    else
        @expression(model, expression_line_losses_nodal[n=SETS.NODES, t=SETS.TIME], 0 )
    end
    #-- Load serving
    # Nodal electricity exports
    @constraint(model, cons_load_serving_elec_export[n=SETS.NODES, t=SETS.TIME],
        p_g_node[n,t] + p_chp_node[n,t] + p_re[n,t] + p_hvdc_node[n,t] + p_s[n,t] + p_bio_node[n,t]
        - p_pth_node[n,t] - p_ptg_node[n,t] - p_d[n,t] - p_hp_node[n,t] - p_emob_node[n,t] - expression_line_losses_nodal[n,t] == p_ex[n,t] + p_slack_pos[n,t] - p_slack_neg[n,t])
    # Load serving
    @constraint(model, cons_load_serving[t=SETS.TIME], sum(p_ex[n,t] for n in SETS.NODES) == 0)

    #-- Load
    # Assignment
    @constraint(model, cons_load_assignment[n=SETS.NODES, t=SETS.TIME], p_d[n,t] == ASSIGNMENTS.A_L[n,:]' * LOADS.TIMESERIES[t, :])

    #-- Load flow
    # line loading calculated but not constrained for debugging informaton #if PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL+RD" # consider line restrictions if nodal setup, when zonal+RD, then do NOT use line restrictions
    # PTDF definition
    # @constraint(model, cons_line_power_flow[m=SETS.BRANCHES, t=SETS.TIME], p_branch[m,t] == BRANCHES.PTDF[m,:]' * p_ex[:,t])
    # Voltage-angle definition
    @constraint(model, cons_line_power_flow_branches[m=SETS.BRANCHES, t=SETS.TIME], p_branch[m,t] == BRANCHES.DELTA_LF_MATRIX_LINE[m,:]' * delta_lf[:,t])
    @constraint(model, cons_line_power_flow_nodes[n=SETS.NODES, t=SETS.TIME], p_ex[n,t] == BRANCHES.DELTA_LF_MATRIX_NODE[n,:]' * delta_lf[:,t])
    # Slack node
    @constraint(model, cons_slack_node[t=SETS.TIME], delta_lf[PARAMETER_SETTINGS.SLACK_NODE, t] == 0 )

    if PARAMETER_SETTINGS.MODEL_TYPE == "CESM" # consider line restrictions if CESM (nodal) setup, when zonal+RD, do NOT use line restrictions
        # determine which lines  to consider
        BOOLEAN_CONSIDER_BRANCH = (BRANCHES.PMAX .<= 9999999)
        #println(BOOLEAN_CONSIDER_BRANCH)
        println("relevant branches for constraints: " * string(sum(BOOLEAN_CONSIDER_BRANCH)) * " out of total amount of: " * string(length(BOOLEAN_CONSIDER_BRANCH)))
        # Max. line capacity (positive)
        @constraint(model, cons_line_max_capacity[m=SETS.BRANCHES, t=SETS.TIME], BOOLEAN_CONSIDER_BRANCH[m] * p_branch_abs[m,t] <= BRANCHES.PMAX[m])
        # # Max. line capacity 2 (negative)
        # @constraint(model, cons_line_max_capacity_2[m=SETS.BRANCHES, t=SETS.TIME], BOOLEAN_CONSIDER_BRANCH[m] * -p_branch[m,t] <= BRANCHES.PMAX[m])
    end

    #-- HVDC

    if PARAMETER_SETTINGS.HVDC_LOSSES == "ACTIVE"
        # HVDC node assignment
        @constraint(model, cons_hvdc_node[n=SETS.NODES, t=SETS.TIME], p_hvdc_node[n,t] == ASSIGNMENTS.A_HVDC[n,:]' * p_hvdc_line[:,t]
            + ASSIGNMENTS.A_HVDC_LOSS[n,:]' * p_hvdc_loss[:,t])
        # Absolute value of HVDC Loading 1
        @constraint(model, cons_hvdc_abs_1[h=SETS.HVDCS, t=SETS.TIME], p_hvdc_line[h,t] <= p_hvdc_abs[h,t])
        # Absolute value of HVDC Loading 2
        @constraint(model, cons_hvdc_abs_2[h=SETS.HVDCS, t=SETS.TIME], -p_hvdc_line[h,t] <= p_hvdc_abs[h,t])
        # Max. line capacity, bounding absolute value of HVDC
        @constraint(model, cons_hvdc_abs_3[h=SETS.HVDCS, t=SETS.TIME], p_hvdc_abs[h,t] <= HVDCS.P_MAX[h])
        # HVDC losses
        @constraint(model, cons_hvdc_losses[h=SETS.HVDCS, t=SETS.TIME], p_hvdc_loss[h,t] == HVDCS.l1_proportional_loss[h] * p_hvdc_abs[h,t] + HVDCS.l0_const_loss[h])
    else
        # Max. HVDC capacity 1 (positive)
        @constraint(model, cons_hvdc_max_capacity_1[h=SETS.HVDCS, t=SETS.TIME], p_hvdc_line[h,t] <= HVDCS.P_MAX[h])
        # Max. line capacity 2  (negative)
        @constraint(model, cons_hvdc_max_capacity_2[h=SETS.HVDCS, t=SETS.TIME], -p_hvdc_line[h,t] <= HVDCS.P_MAX[h])
        # HVDC node assignment
        @constraint(model, cons_hvdc_node[n=SETS.NODES, t=SETS.TIME], p_hvdc_node[n,t] == ASSIGNMENTS.A_HVDC[n,:]' * p_hvdc_line[:,t])
    end



    #-- Conventional power plants
    # Assignment
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_G, p_g_unit, "==", p_g_node, "cons_conventional_assignment");
        #@constraint(model, cons_conventional_assignment[n=SETS.NODES, t=SETS.TIME], p_g_node[n,t] == ASSIGNMENTS.A_G[n,:]' * p_g_unit[:,t])
    # Max. generation capacity
    @constraint(model, cons_conventional_max_capacity[k=SETS.GENS, t=SETS.TIME], p_g_unit[k,t] <= p_g_unit_online[k,t])
    # Max. online capacity
    if ((PARAMETER_SETTINGS.MODEL_TYPE == "INVEST"))
        @constraint(model, cons_conventional_invest[k=SETS.GENS, t=SETS.TIME], p_g_unit_online[k,t] <= p_g_unit_invest[k])
    elseif (PARAMETER_SETTINGS.MODEL_TYPE == "CESM") | (PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL") | (PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL+RD") | (PARAMETER_SETTINGS.MODEL_TYPE == "BENDERS INVEST")
        @constraint(model, cons_conventional_invest[k=SETS.GENS, t=SETS.TIME], p_g_unit_online[k,t] <= CONVENTIONALS.PMAX[k])
    end
    
    # Min. generation capacity
    @constraint(model, cons_conventional_min_capacity[k=SETS.GENS, t=SETS.TIME], CONVENTIONALS.R_G[k] * p_g_unit_online[k,t] <= p_g_unit[k,t] )
    # Costs
    @constraint(model, cons_conventional_costs[k=SETS.GENS, t=SETS.TIME], c_g_op[k,t] == CONVENTIONALS.FUEL_PRICE[k] *
        (CONVENTIONALS.H_OP[k] * p_g_unit[k,t] + (CONVENTIONALS.H_MIN[k] - CONVENTIONALS.H_OP[k]) * CONVENTIONALS.R_G[k] * p_g_unit_online[k,t]))
    # Start-up capacity
    @constraint(model, cons_conventional_startup_capacity[k=SETS.GENS, t=(SETS.TIME[begin]+1):SETS.TIME[end]],
        p_g_unit_online[k,t] - p_g_unit_online[k,t-1] <= p_g_unit_startup[k,t] )
    # @constraint(model, cons_conventional_startup_capacity[k=SETS.CHPS, t=SETS.TIME, p_g_unit_online[k,t] - p_g_unit_online[k,t-1] <= p_g_unit_startup[k,t] )
    # Initial operational state
    @constraint(model, cons_initial_conventional_operational_state[k=SETS.GENS, t=SETS.TIME[begin]],
        p_g_unit_online[k,t]*(1-IS_FIRST_PERIOD) - CONV_ONLINE_LASTPERIOD[k]*(1-IS_FIRST_PERIOD) <= p_g_unit_startup[k,t] )


    #-- Renewables (PV + Wind)
    # Assignment and max. generation capacity
    @constraint(model, cons_renewables_max_capacity[n=SETS.NODES, t=SETS.TIME], p_re[n,t] == p_pv_wind_nodal[n,t] + p_hy_ror_nodal[n,t])
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_REN, p_pv_wind_after_curtailment, "==", p_pv_wind_nodal, "cons_renewables_pv_wind_max_capacity"); # RENEWABLES.TIMESERIES intentionally transposed
    # Max. capacity (and curtailment)
    @constraint(model, cons_pv_wind_infeed[k=SETS.RENS, t=SETS.TIME], p_pv_wind_after_curtailment[k,t] <= RENEWABLES.TIMESERIES[t,k]); # RENEWABLES.TIMESERIES intentionally transposed

    #-- Hydro run-of-river
    # Nodal assignment
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_HY_ROR, p_hy_ror_unit, "==", p_hy_ror_nodal, "cons_renewables_pv_wind_max_capacity"); # RENEWABLES.TIMESERIES intentionally transposed
    # Max. capacity (and curtailment)
    @constraint(model, cons_hy_ror_max[k=SETS.HY_ROR, t=SETS.TIME], p_hy_ror_unit[k,t] <= HY_ROR.TIMESERIES[t,k]) # unnecessary, as curtailment is prevalent


    #-- Biomass/gas
    # Nodal assignment
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_BIO, p_bio_unit, "==", p_bio_node, "cons_bio_assignment");
        #@constraint(model, cons_bio_assignment[n=SETS.NODES, t=SETS.TIME], p_bio_node[n,t] == ASSIGNMENTS.A_BIO[n,:]' * p_bio_unit[:,t])
    # Max. capacity
    @constraint(model, cons_bio_max[k=SETS.BIOS, t=SETS.TIME], p_bio_unit[k,t] <= BIOS.TIMESERIES[t,k])
    # Costs
    @constraint(model, cons_bio_costs[k=SETS.BIOS, t=SETS.TIME], c_bio_op[k,t] == p_bio_unit[k,t] * BIOS.PRICE[k])



    #-- Power to gas  
    if (PARAMETER_SETTINGS.MODEL_TYPE == "INVEST")
        @constraint(model, cons_ptg_invest[k=SETS.PTG, t=SETS.TIME], p_ptg_unit[k,t] <= p_ptg_unit_invest[k])
    elseif (PARAMETER_SETTINGS.MODEL_TYPE == "CESM") | (PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL") | (PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL+RD") | (PARAMETER_SETTINGS.MODEL_TYPE == "BENDERS INVEST")
        @constraint(model, cons_ptg_invest[k=SETS.PTG, t=SETS.TIME], p_ptg_unit[k,t] <= PTG.PMAX[k])
    end
    
    # Assignment to heat regions
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_PTG, p_ptg_unit, "==", p_ptg_node, "cons_ptg_assignment");
    # Hydrogen production
    @constraint(model, cons_ptg_hydrogen_production[k=SETS.PTG, t=SETS.TIME], p_ptg_unit[k,t] * PTG.ETA_P_MAX[k] == p_ptg_unit_hydrogen[k,t] )
    # Load serving
    @constraint(model, cons_ptg_load_serving[t=SETS.TIME], sum(p_ptg_unit_hydrogen[k,t] for k in SETS.PTG) + p_ptg_slack[t] == sum(PTG.PMAX[k] for k in SETS.PTG) )

    # Max. hydrogen capacity
    # @constraint(model, cons_ptg_online_capacity_hydrogen[k=SETS.PTG, t=SETS.TIME], p_ptg_unit_hydrogen[k,t] <= p_ptg_unit_hydrogen_online[k,t] )
    # Hydrogen production
    # @constraint(model, cons_ptg_max_capacity_hydrogen[k=SETS.PTG, t=SETS.TIME], p_ptg_unit_hydrogen_online[k,t] <= PTG.PMAX[k] * PTG.ETA_P_MAX[k])
    # Hydrogen production
    # @constraint(model, cons_ptg_hydrogen_production[k=SETS.GENS, t=SETS.TIME], p_ptg_unit[k,t] == PTG.H_OP[k] * p_ptg_unit_hydrogen[k,t] +
    #     (PTG.H_MIN[k] - PTG.H_OP[k]) * PTG.R_G[k] * p_ptg_unit_hydrogen_online[k,t])



    #-- Combined heat power plants (CHP)
    # Assignment (electricity)
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_CHP_EL, p_chp_unit, "==", p_chp_node, "cons_chp_assignment");
    # Assignment (heat)
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_CHP_H, q_chp_unit, "==", q_chp_region, "cons_chp_assignment_heat");
    # Max. power capacity online
    @constraint(model, cons_chp_max_online_capacity[k=SETS.CHPS, t=SETS.TIME], p_chp_unit_online[k,t] <= COMBINEDHEATPOWERS.PMAX[k])
    # Max. power production
    @constraint(model, cons_chp_max_capacity_bp[k=SETS.CHPS_BP, t=SETS.TIME], p_chp_unit[k,t] <= p_chp_unit_online[k,t])
    @constraint(model, cons_chp_max_capacity_ec[k=SETS.CHPS_EC, t=SETS.TIME], p_chp_unit[k,t] <= p_chp_unit_online[k,t] - COMBINEDHEATPOWERS.CV[k] * q_chp_unit[k,t])
    # Max. heat production
    @constraint(model, cons_chp_max_heat_capacity[k=SETS.CHPS, t=SETS.TIME], q_chp_unit[k,t] <= COMBINEDHEATPOWERS.QMAX[k])
    # Power to heat ratio
    @constraint(model, cons_chp_bp_heatcoupling_bp[k=SETS.CHPS_BP, t=SETS.TIME], COMBINEDHEATPOWERS.CB[k] * q_chp_unit[k,t] == p_chp_unit[k,t])
    @constraint(model, cons_chp_ec_heatcoupling_ec[k=SETS.CHPS_EC, t=SETS.TIME], COMBINEDHEATPOWERS.CB[k] * q_chp_unit[k,t] <= p_chp_unit[k,t])
    # Min. generation capacity
    @constraint(model, cons_chp_min_capacity_bp[k=SETS.CHPS_BP, t=SETS.TIME], COMBINEDHEATPOWERS.R_G[k] * p_chp_unit_online[k,t] <= p_chp_unit[k,t])
    @constraint(model, cons_chp_min_capacity_ec[k=SETS.CHPS_EC, t=SETS.TIME], COMBINEDHEATPOWERS.R_G[k] * p_chp_unit_online[k,t] - COMBINEDHEATPOWERS.CV[k] * q_chp_unit[k,t] <= p_chp_unit[k,t])
    # Costs
    @constraint(model, cons_chp_costs[k=SETS.CHPS, t=SETS.TIME], c_chp_op[k,t] == COMBINEDHEATPOWERS.FUEL_PRICE[k] *
        (COMBINEDHEATPOWERS.H_OP[k] * p_chp_unit[k,t] + (COMBINEDHEATPOWERS.H_MIN[k] - COMBINEDHEATPOWERS.H_OP[k]) *
        COMBINEDHEATPOWERS.R_G[k] * p_chp_unit_online[k,t] + COMBINEDHEATPOWERS.H_OP[k] * COMBINEDHEATPOWERS.CB[k] * q_chp_unit[k,t]))
    # Start-up capacity
    @constraint(model, cons_chp_startup_capacity[k=SETS.CHPS, t=(SETS.TIME[begin]+1):SETS.TIME[end]],
        p_chp_unit_online[k,t] - p_chp_unit_online[k,t-1] <= p_chp_unit_startup[k,t] )
    # @constraint(model, cons_chp_startup_capacity[k=SETS.CHPS, t=SETS.TIME, p_chp_unit_online[k,t] - p_chp_unit_online[k,t-1] <= p_chp_unit_startup[k,t] )
    # Initial operational state
    @constraint(model, cons_initial_chp_operational_state[k=SETS.CHPS, t=SETS.TIME[begin]],
        p_chp_unit_online[k,t]*(1-IS_FIRST_PERIOD) - CHP_ONLINE_LASTPERIOD[k]*(1-IS_FIRST_PERIOD) <= p_chp_unit_startup[k,t] )

    # -- Power to heat
    # Max. power capacity online
    @constraint(model, cons_pth_max_capacity[k=SETS.PTH, t=SETS.TIME], p_pth_unit[k,t] <= PTH.PMAX[k])
    # Assignment to heat regions
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_PTH_H, p_pth_unit, "==", p_pth_region, "cons_pth_heat_assignment");
        #@constraint(model, cons_pth_heat_assignment[r=SETS.HEATREGIONS, t=SETS.TIME], p_pth_region[r,t] == ASSIGNMENTS.A_PTH_H[r,:]' * p_pth_unit[:,t])
    # Assignment to elec_nodes
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_PTH_E, p_pth_unit, "==", p_pth_node, "cons_pth_elec_assignment");
        #@constraint(model, cons_pth_elec_assignment[n=SETS.NODES, t=SETS.TIME], p_pth_node[n,t] == ASSIGNMENTS.A_PTH_E[n,:]' * p_pth_unit[:,t])

    # -- Heat only
    # Max. power capacity online
    @constraint(model, cons_heat_only_max_capacity[k=SETS.HEAT_ONLY, t=SETS.TIME], p_heat_only_unit[k,t] <= HEAT_ONLY.PMAX[k])
    # Assignment
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_HEAT_ONLY, p_heat_only_unit, "==", p_heat_only_region, "cons_heat_only_assignment");
        #@constraint(model, cons_heat_only_assignment[r=SETS.HEATREGIONS, t=SETS.TIME], p_heat_only_region[r,t] == ASSIGNMENTS.A_HEAT_ONLY[r,:]' * p_heat_only_unit[:,t])

    # -- Constant heat load
    # Assignment
        #@constraint(model, cons_const_heat_load_assignment[r=SETS.HEATREGIONS, t=SETS.TIME], p_const_heat_load_region[r,t] == ASSIGNMENTS.A_CONST_HEAT_LOAD[r,:]' * CONST_HEAT_LOAD.Q_MAX)
    @constraint(model, cons_const_heat_load_assignment[r=SETS.HEATREGIONS, t=SETS.TIME], p_const_heat_load_region[r,t] == ASSIGNMENTS.A_CONST_HEAT_LOAD[r,:]' * CONST_HEAT_LOAD.Q_MAX)

    # -- Heat demand
    # Assignment
        
    HEAT_REGIONS_RELEVANT =  SETS.HEATREGIONS[in.(collect(SETS.HEATREGIONS), (collect(1:size(ASSIGNMENTS.A_HEAT_DEMAND,1)),))]
    HEAT_REGIONS_NOT_RELEVANT =  SETS.HEATREGIONS[.!(in.(collect(SETS.HEATREGIONS), (collect(1:size(ASSIGNMENTS.A_HEAT_DEMAND,1)),)))]  
    @constraint(model, cons_heat_demand_assignment[r=HEAT_REGIONS_RELEVANT, t=SETS.TIME], q_heat_demand[r,t] == ASSIGNMENTS.A_HEAT_DEMAND[r,:]' * HEAT_DEMAND.TIMESERIES[t,:])
    @constraint(model, cons_heat_demand_assignment_2[r=HEAT_REGIONS_NOT_RELEVANT, t=SETS.TIME], q_heat_demand[r,t] == 0)

    #@constraint(model, cons_heat_demand_assignment[r=SETS.HEAT_REGIONS, t=SETS.TIME], q_heat_demand[r,t] == ASSIGNMENTS.A_HEAT_DEMAND[r,:]' * HEAT_DEMAND.TIMESERIES[t,:])
    #zum Vergleich: @constraint(model, cons_load_assignment[n=SETS.NODES, t=SETS.TIME], p_d[n,t] == ASSIGNMENTS.A_L[n,:]' * LOADS.TIMESERIES[t, :])


    #-- Heat regions
    # Load serving
    @constraint(model, cons_heatload_serving[r=SETS.HEATREGIONS, t=SETS.TIME], q_chp_region[r,t] + q_chp_region_slack[r,t] +
        p_pth_region[r,t] + p_heat_only_region[r,t] == q_heat_demand[r,t] + p_const_heat_load_region[r,t])





    #-- Heat pumps
    # If timeseries model is selected for flexibilities
    if PARAMETER_SETTINGS.FLEX_MARKET.HP == "TS"
        # Assignment: heat pump to electricity node
        #add_constraint_sparse_2D(model, ASSIGNMENTS.A_HP, HEATPUMPS.ELEC_TIMESERIES', "==", p_hp_node, "cons_heatpump_assignment");
        # @constraint(model, cons_heatpump_assignment[n=SETS.NODES, t=SETS.TIME], p_hp_node[n,t] == ASSIGNMENTS.A_HP[n,:]' * HEATPUMPS.ELEC_TIMESERIES[t,:])
        NODES_RELEVANT =  SETS.NODES[in.(collect(SETS.NODES), (collect(1:size(ASSIGNMENTS.A_HP,1)),))]
        NODES_NOT_RELEVANT =  SETS.NODES[.!(in.(collect(SETS.NODES), (collect(1:size(ASSIGNMENTS.A_HP,1)),)))]  
        @constraint(model, cons_heatpump_assignment[n=NODES_RELEVANT, t=SETS.TIME], p_hp_node[n,t] == ASSIGNMENTS.A_HP[n,:]' * HEATPUMPS.ELEC_TIMESERIES[t,:])
        @constraint(model, cons_heatpump_assignment_zero[n=NODES_NOT_RELEVANT, t=SETS.TIME], p_hp_node[n,t] == 0)
        @expression(model, expression_hp_final_filling_level, 0)
     
    # If detailed model is selected for flexibilities
    elseif PARAMETER_SETTINGS.FLEX_MARKET.HP == "Flex"
    #- Detailed model for heat pumps
        @constraint(model, cons_hp_and_he[s=SETS.HPS, t=SETS.TIME], p_hp_and_he_unit[s,t] == p_hp_unit[s,t] + p_hp_he_unit[s,t] )
        # Assignment: heat pump to electricity node
        add_constraint_sparse_2D(model, ASSIGNMENTS.A_HP, p_hp_and_he_unit, "==", p_hp_node, "cons_heatpump_elecnode_assignment");
        # @constraint(model, cons_heatpump_elecnode_assignment[n=SETS.NODES, t=SETS.TIME], p_hp_node[n,t] == ASSIGNMENTS.A_HP[n,:]' * (p_hp_unit[:,t] + p_hp_he_unit[:,t]))
        # Assignment: heat demand timeseries to heat pump
        @constraint(model, cons_heatpump_heatdemand_assignment[s=SETS.HPS, t=SETS.TIME], q_hp_d[s,t] == HEATPUMPS.HEAT_TIMESERIES[t,s] - hp_heat_slack[s,t])
        # Max. power heat pump
        @constraint(model, cons_heatpump_max_power[s=SETS.HPS, t=SETS.TIME], p_hp_unit[s,t] <= HEATPUMPS.P_MAX_HP[s])
        # Max. power heat element
        @constraint(model, cons_heatpump_heatelement_max_power[s=SETS.HPS, t=SETS.TIME], p_hp_he_unit[s,t] <= HEATPUMPS.P_MAX_HE[s])
        # Filling level
        @constraint(model, cons_heatstorage_filling_level[s=SETS.HPS, t=SETS.TIME],
            q_hp_hs[s,t] == q_hp_hs[s,t-1] - q_hp_d[s,t] +
            ((p_hp_unit[s,t] * HEATPUMPS.COP_TIMESERIES[t,s]) + p_hp_he_unit[s,t]) -
            (q_hp_hs[s,t] + q_hp_hs[s,t-1]) / 2 * HEATPUMPS.NU_SELF_DIS[s]
        )
        # Initial filling level
        @constraint(model, cons_heatpump_heatstorage_initial_filling_level[s=SETS.HPS, t=[(SETS.TIME[begin]-1)]],
            q_hp_hs[s,t]*(1-IS_FIRST_PERIOD) == HP_Q_0[s]*(1-IS_FIRST_PERIOD))
        # Max. filling level
        @constraint(model, cons_heatpump_heatstorage_max_capacity[s=SETS.HPS, t=(SETS.TIME[begin]-1):SETS.TIME[end]], q_hp_hs[s,t] <= HEATPUMPS.Q_MAX[s])
        # Final filling level value
        @expression(model, expression_hp_final_filling_level, sum(HEATPUMPS.PENALTY_LAST_TIMESTEP * q_hp_hs[s,end] for s in SETS.HPS))
    end

    #-- Storage units
    # Assignment
    @constraint(model, cons_storages_assignment_chdis[k=SETS.STORS, t=SETS.TIME], p_dis_minus_ch[k,t] == p_dis[k,t] - p_ch[k,t])
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_STOR, p_dis_minus_ch, "==", p_s, "cons_storages_assignment");
        #@constraint(model, cons_storages_assignment[n=SETS.NODES, t=SETS.TIME], p_s[n,t] == ASSIGNMENTS.A_STOR[n,:]' * (p_dis[:,t] - p_ch[:,t]))
    # Max. filling level
    @constraint(model, cons_storages_max_capacity[s=SETS.STORS, t=(SETS.TIME[begin]-1):SETS.TIME[end]], e_s[s,t] <= STORAGES.EMAX[s])
    # Max. charging power
    @constraint(model, cons_storages_max_charging[s=SETS.STORS, t=SETS.TIME], p_ch[s,t] <= STORAGES.P_MAX_CH[s])
    # Max. charging power
    @constraint(model, cons_storages_max_discharging[s=SETS.STORS, t=SETS.TIME], p_dis[s,t] <= STORAGES.P_MAX_DIS[s])
    # Initial filling level
    @constraint(model, cons_initial_filling_level[s=SETS.STORS, t=[(SETS.TIME[begin]-1)]], e_s[s,t]*(1-IS_FIRST_PERIOD) == E_0[s]*(1-IS_FIRST_PERIOD))
    # Filling level
    @constraint(model, cons_storages_filling_level[s=SETS.STORS, t=SETS.TIME], e_s[s,t] == e_s[s,t-1] + (p_ch[s,t]* STORAGES.EFF_CH[s] - p_dis[s,t]/STORAGES.EFF_DIS[s])
    - (e_s[s,t] + e_s[s,t-1])/2 * STORAGES.SELF_DIS_RATIO[s])
    # Circle constraint
    @constraint(model, cons_storages_circle[s=SETS.STORS], e_s[s,SETS.TIME[end]]*IS_FIRST_PERIOD >= e_s[s,(SETS.TIME[begin]-1)]*IS_FIRST_PERIOD)
    

    #-- Electric Mobility
    # Assignment
    @constraint(model, cons_emob_assignment[n=SETS.NODES, t=SETS.TIME], p_emob_node[n,t] == p_emob_node_residential[n,t] + p_emob_node_commercial[n,t] + p_emob_node_quick[n,t])
    if PARAMETER_SETTINGS.FLEX_MARKET.EMOB == "Flex" 
        # upper boundary on nodal level
        add_constraint_sparse_2D(model, ASSIGNMENTS.A_EMOB_RESIDENTIAL, EMOB_RESIDENTIAL.PMAX', "==", p_emob_node_residential_max, "cons_emob_residential");
        # nodal power must be below upper boundary
        @constraint(model, cons_emob_node_resiential_val[n=SETS.NODES, t=SETS.TIME], p_emob_node_residential[n,t] <= p_emob_node_residential_max[n,t])
        # emob (dumb) demand for time period considered (based on uncontrolled charging timeseries that is grid friendly)
        @constraint(model, cons_emob_residential_dem[n=SETS.NODES, t=SETS.TIME], ASSIGNMENTS.A_EMOB_RESIDENTIAL[n,:]' * EMOB_RESIDENTIAL.TIMESERIES[t,:] == p_emob_node_residential_dumb_dem[n,t]);
        # dumb sum must match Flex timeseries sum
        @constraint(model, cons_emob_residential_energy_sum[n=SETS.NODES], sum(p_emob_node_residential[n,t] for t in SETS.TIME) == sum(p_emob_node_residential_dumb_dem[n,t] for t in SETS.TIME));
    else # "TS"
        add_constraint_sparse_2D(model, ASSIGNMENTS.A_EMOB_RESIDENTIAL, EMOB_RESIDENTIAL.TIMESERIES', "==", p_emob_node_residential, "cons_emob_residential");
    end
    
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_EMOB_COMMERCIAL, EMOB_COMMERCIAL.TIMESERIES', "==", p_emob_node_commercial, "cons_emob_commercial");
    add_constraint_sparse_2D(model, ASSIGNMENTS.A_EMOB_QUICK, EMOB_QUICK.TIMESERIES', "==", p_emob_node_quick, "cons_emob_quick");

  
    # Investment cost
    if ((PARAMETER_SETTINGS.MODEL_TYPE == "INVEST"))
        @expression(model, expression_invest, 
            c_invest_gas * sum(p_g_unit_invest[k] for k in SETS.GENS)
            + c_invest_ptg * sum(p_ptg_unit_invest[k] for k in SETS.PTG))
    elseif ((PARAMETER_SETTINGS.MODEL_TYPE == "CESM") | (PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL") | (PARAMETER_SETTINGS.MODEL_TYPE == "ZONAL+RD")  | (PARAMETER_SETTINGS.MODEL_TYPE == "BENDERS INVEST"))
        @expression(model, expression_invest, 0)
    end
      
    # Objective function
    @objective(model, Min,
        # Fuel prices
        sum(c_g_op[k,t] for k in SETS.GENS, t in SETS.TIME) # Costs of conventionals
        + sum(c_chp_op[k,t] for k in SETS.CHPS, t in SETS.TIME) # Costs of CHP
        + sum(p_bio_unit[k,t] * BIOS.PRICE[k] for k in SETS.BIOS, t in SETS.TIME) # Costs of Biomass/gas units
        + sum(p_heat_only_unit[k,t] * HEAT_ONLY.FUEL_PRICE[k] for k in SETS.HEAT_ONLY, t in SETS.TIME) # Costs of heat only units
        # Start-up costs
        + sum(p_g_unit_startup[k,t]*CONVENTIONALS.START_UP_COSTS[k] for k in SETS.GENS, t in (SETS.TIME[begin]+1):SETS.TIME[end])
        + sum(p_chp_unit_startup[k,t]*COMBINEDHEATPOWERS.START_UP_COSTS[k] for k in SETS.CHPS, t in (SETS.TIME[begin]+1):SETS.TIME[end])
        # Power to gas slack
        + sum(p_ptg_slack[t] * PTG.HYDROGEN_PRICE[1] for t in SETS.TIME)
        
        # Expression for investment decisions
        + expression_invest
        
        # Slacks
        + 99999999 * sum(q_chp_region_slack[r,t] for r in SETS.HEATREGIONS, t in SETS.TIME)
        + 99999999 * sum(p_slack_pos[n,t] for n in SETS.NODES, t in SETS.TIME)
        + 99999999 * sum(p_slack_neg[n,t] for n in SETS.NODES, t in SETS.TIME)
        + 99999999 * sum(hp_heat_slack[s,t] for s in SETS.HPS, t in SETS.TIME)

        # Assign value PENALTY_LAST_TIMESTEP to filling level of last timestep
        - expression_hp_final_filling_level
        - sum(STORAGES.PENALTY_LAST_TIMESTEP * e_s[s,end] for s in SETS.STORS)
    )

    println("---\n", now(), " : market model built has ", num_variables(model), " variables")

    return model

end

function optimize_market_model(model, PARAMETER_SETTINGS, SETS)

    # Optimize
    try
        println("---\n", now(), " : market model: optimize...")
        optimize!(model)
    catch
        println("---\n", now(),  " : market model with presolve = ", string(PARAMETER_SETTINGS.PRESOLVE), " failed. Reverse presolve setting.")
        set_optimizer_attribute(model, "Presolve", abs(PARAMETER_SETTINGS.PRESOLVE-1) )
        optimize!(model)
    end

    #--- Save results in struct
    RESULTS = data_containers.single_model_result(PARAMETER_SETTINGS, SETS, model)


    #--- Return values
    return RESULTS
end
