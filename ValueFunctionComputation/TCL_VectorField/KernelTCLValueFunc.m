classdef KernelTCLValueFunc < TCLValueFunc
    %UNTITLED7 Summary of this class goes here
    %   Detailed explanation goes here

    properties   
        ambiguity_type % Type of ambiguity set
        center_ball % Radius of the mean describing the ambiguity set
        radius_ball % Radius of the variance describing the ambiguity set

        kernel_func
        kernel_parameter

        type_value_func_computation
    end

    properties (Access = private)
        
        % These parameters are only used for the kernel mean embedding
        chol_fac 
        
        number_of_points_KME
        points_to_estimate_RHKS_norm
        chol_gram_matrix_RKHS_norm

        regulariser_param
        eta_param
        data_KME

    end

    methods
        function obj = KernelTCLValueFunc(number_of_points,time_horizon,...
                            type_vector_field,radius_ball,kernel_func,...
                                kernel_parameter,type_value_func_computation,param)
            
            obj = obj@TCLValueFunc(number_of_points,time_horizon,type_vector_field,param);

            obj.ambiguity_type = 'KernelAmbiguity';
            
            obj.type_value_func_computation = type_value_func_computation;

            obj.chol_fac = [];
            obj.number_of_points_KME = [];
            obj.regulariser_param = [];
            obj.eta_param = [];
            obj.data_KME = [];
            obj.points_to_estimate_RHKS_norm = [];

            if strcmp(obj.type_value_func_computation,'KME')
                obj.chol_fac = param.chol_fac;
                obj.number_of_points_KME = param.number_of_points_KME;
                obj.regulariser_param = param.regulariser_param;
                obj.eta_param = param.eta_param;
                obj.kernel_func = param.kernel_func;
                obj.kernel_parameter = param.kernel_parameter;
                obj.data_KME = param.data_KME;

                obj.points_to_estimate_RHKS_norm = linspace(obj.param.safe_set(1),obj.param.safe_set(2),number_of_points)';
                obj.chol_gram_matrix_RKHS_norm = GaussianKernel(obj.points_to_estimate_RHKS_norm',...
                                                    obj.points_to_estimate_RHKS_norm',obj.kernel_parameter,...
                                                    obj.regulariser_param,'MatrixFac',[]);

            else
                obj.kernel_func = kernel_func;
                obj.kernel_parameter = kernel_parameter;
            end

            obj.center_ball = [];
            obj.radius_ball = radius_ball;
            
        end

        function out = iterate_value_function(obj,current_state,current_input...
                                                ,next_value_func,state_partition,norm_value_func_RHKS)
            
            % This function performs a pass on the partition of the state
            % space and computes the corresponding value function


            string_current_state = create_x_and_u_string(current_state,[]);

            index_current_state = strcmp(obj.param.grid_no_inputs,string_current_state);

            index_safety = obj.index_safe_set;

            if ~index_safety(index_current_state)
                out = 0; % Safety probability equal to zero if outside the safe set
            else
                out = obj.inner_optimisation(current_state,current_input,...
                                            next_value_func,state_partition,norm_value_func_RHKS);
            end
        end
        
        function out = inner_optimisation(obj,current_state,current_input,...
                                                next_value_func,state_partition,norm_value_func_RHKS)
            % Returns the value function for a given state-action pair (current_state,current_input) using the
            % value function at the next iteration (next_value_func).

            if isempty(obj.type_value_func_computation)
                error(['This function can only be used after initialising ' ...
                    'the field type_value_func_computation to one of these values:' ...
                        '   1. Conservative \n 2. Matrix \n 3. QP \n 4. KME'])
            end
            
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

            obj.center_ball = trans_current_state_input{1};
            grid_x = state_partition.partition.grid_x;

            if strcmp(obj.type_value_func_computation,'KME')
               TCL_kernel_obj = KernelBasedAmbiguity(objective_cost,obj.radius_ball,...
                                                    obj.center_ball,obj.kernel_func,...
                                                         obj.kernel_parameter,obj.param.eta_param,grid_x,...
                                                            obj.type_value_func_computation,...
                                                               obj.chol_fac,norm_value_func_RHKS);
            else
                TCL_kernel_obj = KernelBasedAmbiguity(objective_cost,obj.radius_ball,...
                                                    obj.center_ball,obj.kernel_func,...
                                                         obj.kernel_parameter,grid_x,...
                                                            obj.type_value_func_computation,...
                                                            [],[]);
            end

            

            % Ambiguity parameters
            TCL_kernel_obj.number_of_samples_KME = obj.param.outer_loop_info.string_ambiguity{1}.number_of_samples_KME;

            TCL_kernel_obj.current_state = current_state;
            TCL_kernel_obj.current_input = current_input; 
            
            switch obj.type_value_func_computation
                case 'Conservative'
                    TCL_kernel_obj = TCL_kernel_obj.solve_optimisation_conservative(state_partition);
                    out = TCL_kernel_obj.results_optimisation.opt_obj;
                case 'Matrix'
                    TCL_kernel_obj = TCL_kernel_obj.solve_optimisation(state_partition);
                    out = TCL_kernel_obj.results_optimisation.opt_obj;
                case 'QP'
                    TCL_kernel_obj = TCL_kernel_obj.solve_optimisation_QP(state_partition);
                    out = TCL_kernel_obj.results_optimisation.opt_obj;
                case 'KME'
                    TCL_kernel_obj = TCL_kernel_obj.solve_optimisation_KME(state_partition,obj.data_KME);
                    out = TCL_kernel_obj.results_optimisation.opt_obj;
                otherwise
                    not_implemented();
            end
   
        end

        function obj = backward_iteration(obj,state_partition,input_partition,...
                                            outer_loop_info)
            
            if isempty(obj.type_value_func_computation)
                error(['This function can only be used after initialising ' ...
                    'the field type_value_func_computation to one of these values:' ...
                        '   1. Conservative \n 2. Matrix \n 3. QP \n 4. KME'])
            end


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
            total_iterations = double(time_horizon)*(number_of_points)*number_of_inputs;
            print_inner_loop(total_iterations,0,0,obj.ambiguity_type,outer_loop_info);

            for i = time_horizon:-1:1
                next_value_func = value_function(:,i+1); % saving in a temporary variable the 
                                                       % value function of the next step

                norm_value_func_RHKS = compute_norm_RKHS(next_value_func,obj.chol_gram_matrix_RKHS_norm);

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
                                   current_input,next_value_func,state_partition,norm_value_func_RHKS); % getting the new value for 
                                                                                   % the value function
                    end
                    [value_func_temp(j),opt_input_temp(j)] = max(value_func_temp_input);
                end

                final_time = toc(initial_time);
                

                value_function(:,i) = value_func_temp;
                opt_input(:,i) = opt_input_temp;

                % Printing on the screen
                temp_int = double(time_horizon-i+1);
                iterates = remaining_iterations(1,[temp_int,double(time_horizon)],number_of_points*number_of_inputs,[]); % This is the number of iterations 
                                                                                                                 % completes so far. The name of the matlab 
                                                                                                                 % function may be misleading

                print_inner_loop(total_iterations,iterates,final_time,obj.ambiguity_type,outer_loop_info);
            end
            
            
            obj.value_function = value_function;
            obj.opt_input = opt_input;

        end
    end
end

function norm_RKHS = compute_norm_RKHS(value_function,chol_factorization)

beta = chol_factorization\(chol_factorization'\value_function);
norm_RKHS = sqrt(beta'*chol_factorization'*chol_factorization*beta);


end