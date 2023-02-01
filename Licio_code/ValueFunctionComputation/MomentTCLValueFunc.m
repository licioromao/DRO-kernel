classdef MomentTCLValueFunc < TCLValueFunc
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ambiguity_type % Type of ambiguity set
        radius_mean % Radius of the mean describing the ambiguity set
        radius_variance % Radius of the variance describing the ambiguity set
    end

    methods
        function obj = MomentTCLValueFunc(number_of_points,time_horizon,...
                            type_vector_field,radius_mean,radius_variance,param)
            
            obj = obj@TCLValueFunc(number_of_points,time_horizon,type_vector_field,param);

            obj.ambiguity_type = 'MomentAmbiguity';

            obj.radius_mean = radius_mean;
            obj.radius_variance = radius_variance;
        end

        function out = iterate_value_function(obj,current_state,current_input...
                                                ,next_value_func,state_partition)
            
            % This function performs a pass on the partition of the state
            % space and computes the corresponding value function


            string_current_state = create_x_and_u_string(current_state,[]);

            index_current_state = strcmp(obj.param.grid_no_inputs,string_current_state);

            index_safety = obj.index_safe_set;

            if ~index_safety(index_current_state)
                out = 0; % Safety probability equal to zero if outside the safe set
            else
                out = obj.inner_optimisation(current_state,current_input,...
                                            next_value_func,state_partition);
            end
        end
        
        function out = inner_optimisation(obj,current_state,current_input,...
                                                next_value_func,state_partition)
            % Returns the value function for a given state-action pair (current_state,current_input) using the
            % value function at the next iteration (next_value_func).
            
            transition_prob = obj.param.transition_prob; % vector with the transition probability matrix
            grid_no_inputs = obj.param.grid_no_inputs;
            
            if length(next_value_func) ~= length(grid_no_inputs)
                error('The dimension of the thrid argument is incorrect.') % outputs an error if 
                                                                           % L is inconsistent with the size 
                                                                           % of the input Z
            end
            
            string_state_input = create_x_and_u_string(current_state,current_input);
            
            if ~transition_prob.isKey({string_state_input}) % outputs an error is there is 
                                                            % no element in transition_prob with the 
                                                            % label (string_currente_state,string_current_input)
                error('The input-action pair is not a member of the transition probability')
            end          
            
            objective_cost = next_value_func;
            trans_current_state_input = transition_prob.values({string_state_input});
            temp_partition = state_partition.get_values.partition;


            [support_set_distribution,mean_center,variance_center] = compute_support_set_mean_variance...
                            (trans_current_state_input{1},temp_partition.grid_x,'WithSigma');

            TCL_moment_obj = MomentBasedAmbiguity(objective_cost,variance_center,mean_center,...
                            obj.radius_variance,obj.radius_mean,support_set_distribution);

            opt_results= TCL_moment_obj.solve_optimisation;

            if (opt_results.results_optimisation.solver_status.problem == 0) ...
                    || (opt_results.results_optimisation.solver_status.problem == 4) ...
                        || (opt_results.results_optimisation.solver_status.problem == -1)
                out = opt_results.results_optimisation.optimal_obj;
            else
                error('There is a problem when solving the optimization problem')
            end

        end

        function obj = backward_iteration(obj,state_partition,input_partition,...
                                            outer_loop_info)


            % Testing the value of N
            if isempty(obj.time_horizon)
                error('Please, initialize the field time_horizon before calling this function');
            elseif obj.time_horizon < 0 || ~isinteger(obj.time_horizon)
                error('The field time_horizon must be a positive integer (int8, int16, etc...)');
            end
            
            % Setting up the optimisation problem
            number_of_points = size(state_partition.partition.grid_x,1); % number of points of the value function
            number_of_inputs = size(input_partition,1);


            if isempty(obj.index_safe_set)
                error(['This function cannot run if index_safe_set is empty. Please try using the ' ...
                                            'method get_index_safety()']);
            end
            
            time_horizon = obj.time_horizon;

            value_function = zeros(number_of_points,time_horizon + 1);
            opt_input = zeros(number_of_points,time_horizon + 1);
            value_function(obj.index_safe_set,time_horizon+1) = 1; % initializing the value function on the safe set

            grid_x = state_partition.get_values.partition.grid_x;

            % Creating a progress bar of the value function computation
            total_iterations = time_horizon*(number_of_points)*number_of_inputs;
            print_inner_loop(total_iterations,0,0,obj.ambiguity_type,outer_loop_info);

            for i = time_horizon:-1:1
                next_value_func = value_function(:,i+1); % saving in a temporary variable the 
                                                       % value function of the next step
                value_func_temp = zeros(number_of_points,1);
                opt_input_temp = zeros(number_of_points,1);

                initial_time = tic;

                for j = 1:number_of_points % iterates over the number of points
                    current_state = grid_x(j,:)'; 

                    value_func_temp_input = zeros(number_of_inputs,1);
                    % Iterating over inputs
                    for u_counter = 1:number_of_inputs
                        current_input = input_partition(u_counter,:)'; 
                        value_func_temp_input(u_counter) = ...
                            obj.iterate_value_function(current_state,...
                                   current_input,next_value_func,state_partition); % getting the new value for 
                                                                                   % the value function
                    end
                    [value_func_temp(j),opt_input_temp(j)] = max(value_func_temp_input);
                end

                final_time = toc(initial_time);
                

                value_function(:,i) = value_func_temp;
                opt_input(:,i) = opt_input_temp;

                % Printing on the screen
                temp_int = double(time_horizon-i+1);
                iterates = remaining_iterations(1,[temp_int,time_horizon],number_of_points*number_of_inputs,[]); % This is the number of iterations 
                                                                                                                 % completes so far. The name of the matlab function may be misleading
                print_inner_loop(total_iterations,iterates,final_time,obj.ambiguity_type,outer_loop_info);
            end

            obj.value_function = value_function;
            obj.opt_input= opt_input;

        end
    end
end


function [support_set_distribution,mean_center,variance_center] = ...
    compute_support_set_mean_variance(trans_prob,grid,type)

% This private method is used by inner approximation to compute
% the center of mean and variance, as well as the support of
% the distribution.

support_set_distribution = grid';
mean_center = support_set_distribution*trans_prob;

number_of_points_variance = 1000;

if strcmp(type,'WithSigma')
    samples = discretesample(trans_prob,number_of_points_variance);
    sum_sigma = zeros(1,1);

    for i =1:number_of_points_variance
        x = support_set_distribution(:,samples(i));
        temp = x - mean_center;
        sum_sigma = temp*temp' + sum_sigma;
    end

    variance_center = sum_sigma/number_of_points_variance;
else
    variance_center = [];
end
end