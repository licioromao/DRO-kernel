classdef NoAmbiguityLTIValueFunc < LTIValueFunc
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ambiguity_type % Type of ambiguity set
    end

    methods
        function obj = NoAmbiguityLTIValueFunc(number_of_points,time_horizon,...
                            type_vector_field,param)
            
            obj = obj@LTIValueFunc(number_of_points,time_horizon,type_vector_field,param);

            obj.ambiguity_type = 'NoAmbiguity';
        end

        function out = iterate_value_function(obj,current_state,current_input...
                                                ,next_value_func)

            % This function performs a pass on the partition of the state
            % space and computes the corresponding value function

            out = obj.inner_optimisation(current_state,current_input,...
                next_value_func);

        end
        
        function out = inner_optimisation(obj,current_state,current_input,...
                                                next_value_func)
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


            out = trans_current_state_input{1}'*objective_cost;


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
            opt_input = zeros(number_of_points,time_horizon+1);
            
            grid_x = state_partition.get_values.partition.grid_x;

            value_function(:,time_horizon+1) = initialise_value_func(grid_x,obj.param.Q); % initializing the value function on the safe set

           
%             % Creating a progress bar of the value function computation
%             total_iterations = double(time_horizon)*(number_of_points)*number_of_inputs;
%             print_inner_loop(total_iterations,0,0,obj.ambiguity_type,outer_loop_info);

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
                                   current_input,next_value_func); % getting the new value for 
                                                                                   % the value function
                    end
                    [value_func_temp(j),opt_input_temp(j)] = max(value_func_temp_input);
                end

                final_time = toc(initial_time);
                

                value_function(:,i) = value_func_temp;
                opt_input(:,i) = opt_input_temp;

%                 % Printing on the screen
%                 temp_int = double(time_horizon-i+1);
%                 iterates = remaining_iterations(1,[temp_int,time_horizon],number_of_points*number_of_inputs,[]); % This is the number of iterations 
%                                                                                                                  % completes so far. The name of the matlab function may be misleading
%                 print_inner_loop(total_iterations,iterates,final_time,obj.ambiguity_type,outer_loop_info);
            end

            obj.value_function = value_function;
            obj.opt_input= opt_input;

        end
    end
end


