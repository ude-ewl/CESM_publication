module basic_functions

using DataFrames, LinearAlgebra, SparseArrays

#--- Function defintions

function create_assignment_matrix(IN_NODE, IN_DF, column_df_name, column_node_name)
    tmp_name = "node_name"
    IN_DF[:,tmp_name] = [findall( x -> x == IN_DF[idx,column_df_name], IN_NODE[:,column_node_name])[1] for idx = 1:size(IN_DF,1)]
    vec_row = Int64[]
    vec_col = Int64[]
    vec_val = Int64[]
    for i in 1:size(IN_DF,1)
        append!(vec_row,IN_DF[i,tmp_name])
        append!(vec_col,i)
        append!(vec_val,1)
    end
    A_MATRIX = sparse(vec_row, vec_col, vec_val,size(IN_NODE,1),size(IN_DF,1))
    return A_MATRIX
end

# Function to append slack node to PTDF matrix
function append_slack_2_ptdf(PTDF_NO_SLACK, SLACK_NODE)
    PTDF_NO_SLACK = DataFrame(PTDF_NO_SLACK, :auto)
    insertcols!(PTDF_NO_SLACK, SLACK_NODE, "x" => zeros(size(PTDF_NO_SLACK,1)))
    return Matrix( PTDF_NO_SLACK)
end

# Function to create assignment matrices
function create_assignment_matrix_timeseries(IN_NODE, IN_DF, IN_DF_TIMESERIES, column_df_name, column_node_name)
    # Check if DataFrames 'IN_DF' and 'IN_DF_TIMESERIES' have the same number of columns
    if size(IN_DF_TIMESERIES, 2)-1 != size(IN_DF, 1)
        error("DataFrames 'IN_DF' and 'IN_DF_TIMESERIES' have different number of columns! Check DataFrames and try again ...")
    end
    # Define column names
    column_node_id = "node_id"
    column_names = names(IN_DF_TIMESERIES)
    # Get index of timeseries column in df
    idx_timeseries_2_df = [findall( x -> x == column_names[idx], IN_DF[:,column_df_name])[1] for idx = 2:size(IN_DF_TIMESERIES,2)]

    df_merge = IN_DF[idx_timeseries_2_df,:]
    # Check if DataFrames 'IN_DF' and 'IN_DF_TIMESERIES' have the same number of columns
    if size(unique(df_merge),1) != size(IN_DF,1)
        error("DataFrame 'IN_DF_TIMESERIES' contains some entities (either loads or renewables) twice! Check DataFrame and try again ...")
    end

    df_merge[:,column_node_id] = [findall( x -> x == df_merge[idx,column_node_name], IN_NODE[:,column_node_name])[1] for idx = 1:size(IN_DF,1)]
    #A_MATRIX = zeros(size(IN_NODE,1), size(IN_DF,1))
    vec_row = Int64[]
    vec_col = Int64[]
    vec_val = Float64[]
    for i in 1:size(IN_DF,1)
        #A_MATRIX[df_merge[i,column_node_id],i] = 1
        append!(vec_row,df_merge[i,column_node_id])
        append!(vec_col,i)
        append!(vec_val,1)
    end
    A_MATRIX = sparse(vec_row, vec_col, vec_val,size(IN_NODE,1),size(IN_DF,1))
    return A_MATRIX
end

function create_PTDF(IN_LINE, IN_NODE, COLUMNS_LINE, COLUMNS_NODE, SLACK_NODE)
    COLUMNS_LINE, COLUMNS_NODE
    # Make (temporary) copy to not mutate/change original data
    temp_IN_LINE = copy(IN_LINE);
    temp_IN_NODE = copy(IN_NODE);
    # Create numerical f_node (from node) and t_node (to node)
    temp_IN_LINE[:,COLUMNS_LINE[2]] = [findall( x -> x == temp_IN_LINE[idx,COLUMNS_LINE[2]], temp_IN_NODE[:,COLUMNS_NODE[1]])[1] for idx = 1:size(IN_LINE,1)]
    temp_IN_LINE[:,COLUMNS_LINE[3]] = [findall( x -> x == temp_IN_LINE[idx,COLUMNS_LINE[3]], temp_IN_NODE[:,COLUMNS_NODE[1]])[1] for idx = 1:size(IN_LINE,1)]

    # Create A Matrix (lines x nodes)
    A = zeros(size(IN_LINE,1), size(IN_NODE,1))
    A_FROM = zeros(size(IN_LINE,1), size(IN_NODE,1))
    A_TO = zeros(size(IN_LINE,1), size(IN_NODE,1))
    for i in 1:size(IN_LINE,1)
        A[i,temp_IN_LINE[i,COLUMNS_LINE[2]]]=1
        A[i,temp_IN_LINE[i,COLUMNS_LINE[3]]]=-1
        A_FROM[i,temp_IN_LINE[i,COLUMNS_LINE[2]]]=1
        A_TO[i,temp_IN_LINE[i,COLUMNS_LINE[3]]]=1
    end

    # Remove slack node from A => create A_
    A_ = A[1:end, 1:end .!=SLACK_NODE]

    # Create B Matrix
    B = Diagonal(temp_IN_LINE[:,COLUMNS_LINE[6]].^(-1))

    # Create PTDF Matrix
    PTDF_NO_SLACK = B * A_ * inv( A_' * B * A_ )

    # Create Matrix for lines
    DELTA_LF_MATRIX_LINE = B*A
    # Create Matrix for nodes
    DELTA_LF_MATRIX_NODE = A'*B*A

    RESISTANCE = temp_IN_LINE[:,COLUMNS_LINE.r]

    # Append slack to PTDF
    return append_slack_2_ptdf(PTDF_NO_SLACK, SLACK_NODE), DELTA_LF_MATRIX_LINE, DELTA_LF_MATRIX_NODE
    return append_slack_2_ptdf(PTDF_NO_SLACK, SLACK_NODE), DELTA_LF_MATRIX_LINE, DELTA_LF_MATRIX_NODE, A_FROM', A_TO', RESISTANCE
end

function create_assignment_matrix_HVDC(IN_LINE, IN_NODE, COLUMNS_NODE, COLUMNS_HVDC)
    # Make (temporary) copy to not mutate/change original data
    temp_IN_LINE = copy(IN_LINE);
    temp_IN_NODE = copy(IN_NODE);
    # Create numerical f_node (from node) and t_node (to node)
    temp_IN_LINE[:,COLUMNS_HVDC.node_id1] = [findall( x -> x == temp_IN_LINE[idx,COLUMNS_HVDC.node_id1], temp_IN_NODE[:,COLUMNS_NODE.node_id])[1] for idx = 1:size(IN_LINE,1)]
    temp_IN_LINE[:,COLUMNS_HVDC.node_id2] = [findall( x -> x == temp_IN_LINE[idx,COLUMNS_HVDC.node_id2], temp_IN_NODE[:,COLUMNS_NODE.node_id])[1] for idx = 1:size(IN_LINE,1)]

    # Create A Matrix (lines x nodes)
    A = zeros(size(IN_LINE,1), size(IN_NODE,1))
    for i in 1:size(IN_LINE,1)
        A[i,temp_IN_LINE[i,COLUMNS_HVDC.node_id1]]=1
        A[i,temp_IN_LINE[i,COLUMNS_HVDC.node_id2]]=-1
    end
    A_loss = zeros(size(IN_LINE,1), size(IN_NODE,1))
    for i in 1:size(IN_LINE,1)
        A_loss[i,temp_IN_LINE[i,COLUMNS_HVDC.node_id1]]=0.5
        A_loss[i,temp_IN_LINE[i,COLUMNS_HVDC.node_id2]]=0.5
    end
    return A', A_loss'
end

#--- Database only

function preprocess_ac_lines_and_transformers(data_elec_network, COLS_TRANSFORMER, COLS_AC_LINE, COLS_NODE, ref_voltage)
    # Names of technical parameters of transformers
    short_circuit_voltage = COLS_TRANSFORMER.short_circuit_voltage
    nominal_apparent_power = COLS_TRANSFORMER.nominal_apparent_power

    # Names of technical parameters of transformers
    reactance_name = COLS_AC_LINE.x
    resistance_name = COLS_AC_LINE.r
    length = COLS_AC_LINE.length
    voltage_level = COLS_AC_LINE.voltage_level
    max_current = COLS_AC_LINE.imax

    # Names of node identifiers
    node_1 = COLS_TRANSFORMER.node_id1
    node_2 = COLS_TRANSFORMER.node_id2
    node_name = COLS_NODE.node_id

    # Internal names
    reactance = "x"
    tmp_reactance = "Reactance"
    resistance = "r"
    tmp_resistance = "Resistance"
    f_node = "node_id1"
    t_node = "node_id2"
    voltage_level_f = "voltage_level_f"
    voltage_level_t = "voltage_level_t"
    max_capacity = "nominal_apparent_power"
    branch_type = "Branch_type"
    transformer = "Transformer"
    ac_line = "Line"
    voltage_level = "voltage_level"

    ## Create PTDF matrix

    in_trafos = data_elec_network.transformers
    in_lines = data_elec_network.ac_lines
    in_nodes = data_elec_network.nodes

    # Number of transformers lines and nodes
    N_lines = size(in_lines,1)
    N_nodes = size(in_nodes,1)
    N_trafos = size(in_trafos,1)

    # Create empty dataframe for output
    branches = DataFrame()

    # Prepare and add transformer information to branches

    # Get reactance (impedance) of transformers: X_k ≈ Z_k = (u_k / S_n) * (U_ref^2 / 100) with Z_k = X_k (approx.)
    in_trafos[:,reactance] = in_trafos[:, short_circuit_voltage] ./ in_trafos[:, nominal_apparent_power] * ref_voltage^2/100
    # Get resistance of transformers: R_k ≈ Z_k * 0.01
    in_trafos[:,resistance] = in_trafos[:, short_circuit_voltage] ./ in_trafos[:, nominal_apparent_power] * ref_voltage^2/100 * 0.01
    # Get node IDs of branches
    in_trafos[:,f_node] = in_trafos[:,node_1] # [findall( x -> x == in_trafos[idx,node_1], in_nodes[:,node_name])[1] for idx = 1:N_trafos]
    in_trafos[:,t_node] = in_trafos[:,node_2] # [findall( x -> x == in_trafos[idx,node_2], in_nodes[:,node_name])[1] for idx = 1:N_trafos]
    # Get voltage levels of branches
    in_trafos[:,voltage_level_f] = -99 *ones(N_trafos)
    in_trafos[:,voltage_level_t] = -99 *ones(N_trafos)
    # Define max. power of transformer in MW
    #print(in_trafos[:, max_capacity])
    in_trafos[:, max_capacity] = in_trafos[:, max_capacity]
    #in_trafos[:, max_capacity] = 9999*ones(N_trafos)
    # Indicate type of branch
    in_trafos[:, branch_type] .= transformer

    append!(branches, DataFrame(x=in_trafos[:,reactance], r=in_trafos[:,resistance], node_id1=in_trafos[:,f_node],
    node_id2=in_trafos[:,t_node], Pmax = in_trafos[:, max_capacity], Voltage_level_f=in_trafos[:,voltage_level_f] , Voltage_level_t=in_trafos[:,voltage_level_t],
    Branch_type = in_trafos[:, branch_type] ))


    # Prepare line information of branches, then append afterwards to branches dataframe

    # Convert reactance of lines to reference voltage level by X_ref = X_l * (U_ref / U_N)^2
    in_lines[:,tmp_reactance] = in_lines[:,reactance_name] .* in_lines[:,length] .*( ref_voltage ./ in_lines[:,voltage_level] ).^2
    # Convert resistance of lines to reference voltage level by R_ref = R_l * (U_ref / U_N)^2
    in_lines[:,tmp_resistance] = in_lines[:,resistance_name] .* in_lines[:,length] .*( ref_voltage ./ in_lines[:,voltage_level] ).^2
    # Get max. power of line in MW
    in_lines[:,max_capacity] = 3^0.5 * in_lines[:,voltage_level] .* in_lines[:,max_current]

    # overwrite low voltage lines with arbitrarily high thermal Limit
    # Get node IDs of lines
    in_lines[:,f_node] = in_lines[:,node_1] # [findall( x -> x == in_lines[idx,node_1], in_nodes[:,node_name])[1] for idx = 1:N_lines]
    in_lines[:,t_node] = in_lines[:,node_2] # [findall( x -> x == in_lines[idx,node_1], in_nodes[:,node_name])[1] for idx = 1:N_lines]
    # Get voltage levels of branches
    in_lines[:,voltage_level_f] = in_lines[:,voltage_level]
    in_lines[:,voltage_level_t] = in_lines[:,voltage_level]

    #in_lines[!,voltage_level_f] = convert.(Float64,in_lines[:,voltage_level_f])
    #in_lines[!,voltage_level_t] = convert.(Float64,in_lines[:,voltage_level_t])

    # Indicate type of branch
    in_lines[:, branch_type] .= ac_line

    # Append lines to branches
    if nrow(branches) > 0 
        append!(branches, DataFrame(x=in_lines[:,tmp_reactance], r=in_lines[:,tmp_resistance], node_id1=in_lines[:,f_node],
            node_id2=in_lines[:,t_node], Pmax = in_lines[:, max_capacity], Voltage_level_f=in_lines[:,voltage_level_f] , Voltage_level_t=in_lines[:,voltage_level_t],
            Branch_type=in_lines[:, branch_type] ))
    else # if model has no transformers (adding empty dataframe "branches" to line dataframe causes errors do to appending missing and non-missing values together)
        branches = DataFrame(x=in_lines[:,tmp_reactance], r=in_lines[:,tmp_resistance], node_id1=in_lines[:,f_node],
        node_id2=in_lines[:,t_node], Pmax = in_lines[:, max_capacity], Voltage_level_f=in_lines[:,voltage_level_f] , Voltage_level_t=in_lines[:,voltage_level_t],
        Branch_type=in_lines[:, branch_type])
    end
    
    return branches, in_trafos, in_lines

end

function get_technical_parameters(IN_DF, COLUMNS_TECHNOLOGY)
    COLS_IN = COLUMNS_TECHNOLOGY
    COLS_OUT = COLUMNS_TECHNOLOGY

    IN_DF[:,COLS_OUT.min_load_factor] = IN_DF[:,COLS_IN.power_elec_min] ./ IN_DF[:,COLS_IN.power_elec_max]
    IN_DF[:,COLS_OUT.heat_rate_elec_min] = ones(size(IN_DF,1)) ./ IN_DF[:,COLS_IN.eta_power_elec_min]
    tmp_h_max = ones(size(IN_DF,1)) ./ IN_DF[:,COLS_IN.eta_power_elec_max]
    tmp_h_min = ones(size(IN_DF,1)) ./ IN_DF[:,COLS_IN.eta_power_elec_min]
    tmp_r = IN_DF[:,COLS_OUT.min_load_factor]

    IN_DF[:,COLS_OUT.heat_rate_marginal] = (tmp_h_max - tmp_h_min .* tmp_r) ./ (ones(size(IN_DF,1)) - tmp_r)
    return IN_DF
end


end
