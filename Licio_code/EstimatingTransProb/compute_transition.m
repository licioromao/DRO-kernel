function param = compute_transition(type_vector_field,number_of_points,...
                                        number_of_MC_simulations)

switch type_vector_field
    case 'TCL'
        %% Model parameters

        R = 2; % Thermal resistence in Celsius/kW
        C = 2; % Thermal capacity ini kW/Celsius
        P = 14; % Range of energy transfer to or from the thermal mass in kW
        eta = 0.7; % Control efficiency
        h = 5/60; % Discretization time in hours

        alpha = exp(-h/(C*R));
        
        grid = [18,23]; % Range of the partition
        Al = 19; % Lower bound on the safe set
        Ah = 22; % Upper bound on the safe set

        mu = 0; sigma = 100*0.25^2; % Estimate on the mean and variance
        W = [-0.5*sqrt(sigma/12),0.5*sqrt(sigma/12)];  % Support of the distribution

        theta = 32; % Environment temperature
        p = 0.05; % success probability
        
        param.grid = grid;
        param.alpha = alpha;
        param.theta = theta;
        param.eta = eta;
        param.R = R;
        param.P = P;
        param.C = C;
        param.h = h;
        param.p = p;
        param.safe_set = [Al,Ah];

        pd = makedist('Normal','mu',mu,'sigma',sigma);
        pd_truncated = truncate(pd,W(1),W(2));

        param.w = pd_truncated;

        param.number_of_points = number_of_points;
        input_partition = generate_input_partition([],'TCL'); % Input partition


        state_partition = OneDimStatePartition(grid,number_of_points,...
                                        param.safe_set,'TCL'); % Partition of the state space

        [grid_with_inputs,grid_no_inputs] = state_partition.create_list(input_partition); % List containing labels for the discrete states of the discretazation

        param.grid_with_inputs = grid_with_inputs;
        param.grid_no_inputs = grid_no_inputs;
        param.size_partition = state_partition.get_size_partition();

        param.number_of_MC_simulations = number_of_MC_simulations;
        
        transition_prob = estimate_transition(state_partition,input_partition,param); 

        % Generating the transition probability
        param.transition_prob = transition_prob;

    otherwise
        NotImplemented();
end




end