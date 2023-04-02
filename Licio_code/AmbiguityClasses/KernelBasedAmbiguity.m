classdef KernelBasedAmbiguity < DistanceBasedAmbiguity
    % This class implements the computation of the Kernel-based ambiguity
    % set problem. It inherits from the class DistanceBasedAmbiguity
    
    % This is very specific of the Gaussian Kernel. Can be generalized later.    
    properties
        param % This variable must be made global at some point
        
        number_of_samples_KME
        norm_value_func_RHKS
        
        current_state
        current_input
        
        kernel_function % This parameter is a handle to a function that defined the Kernel on the state space
        kernel_parameter
        eta_param
        chol_factorization
    end
    
    methods
        function obj = KernelBasedAmbiguity(objective_cost,radius_ball, ... 
                        center_ball,kernel_function,kernel_parameter,eta_param...
                                ,grid_x,type_value_func_computation,chol_fac,norm_value_func_RHKS)
                            

            obj = obj@DistanceBasedAmbiguity(objective_cost,radius_ball,center_ball); % Calling the connstructor of the parent class
                
            obj.param = []; % This should be made a global variable at some point
            
            obj.number_of_samples_KME = [];
           
            obj.current_state = [];
            obj.current_input = [];
            
            % Assigning the kernel parameters to the object
            obj.kernel_parameter = kernel_parameter;
            obj.kernel_function = kernel_function;
            obj.eta_param = eta_param;
            obj.norm_value_func_RHKS = norm_value_func_RHKS;

            if strcmp(type_value_func_computation,'KME')
                obj.chol_factorization = chol_fac;
            else
                obj.chol_factorization = kernel_function(grid_x',grid_x',kernel_parameter,[],'MatrixFac');
            end
            
        end
                       
        function obj = solve_optimisation_conservative(obj,state_partition)
            % This function implements the solution of a kernel ambiguity
            % minimisation where the center of the ambiguity ball is equal
            % to the empirical estimate. 
            % WE MAY WANT TO DELETE THIS FUNCTION AT A LATER VERSION OF THE
            % CODE. THE NUMERICAL RESULTS WITH THIS DOES NOT LOOK
            % PROMISSING. 

            if isempty(obj.kernel_parameter) || isempty(obj.number_of_samples_KME) ... 
                      || isempty(obj.current_state) || isempty(obj.current_input)
                error(['You must first initialize the kernel_parameter and number_of_data ' ...
                    'fiels of the object, as well as current_state and current_input ' ...
                    'using the set_value method']);
            end

            center_ball = obj.center_ball;
            grid_x = state_partition.get_values.partition.grid_x; % CHECK THIS AFTER CHANGING PARTITION OBJECT

            tempKernel = obj.kernel_function(grid_x',grid_x', ... 
                                    obj.kernel_parameter,center_ball,'Normal');
            

            obj.results_optimisation.opt_obj = max(0,center_ball'*obj.objective_cost -... 
                                        (obj.radius_ball)*tempKernel-2*obj.radius_ball);
        end

        function obj = solve_optimisation(obj,state_partition)

            % WE MAY WANT TO DELETE THIS FUNCTION AT A LATER VERSION OF THE
            % CODE. THE NUMERICAL RESULTS WITH THIS DOES NOT LOOK
            % PROMISSING. 
            if isempty(obj.kernel_parameter) || isempty(obj.number_of_samples_KME) ...
                  || isempty(obj.current_state) || isempty(obj.current_input)
                error(['You must first initialize the kernel_parameter and ' ...
                    'number_of_data fiels of the object, as well as ' ...
                    'current_state and current_input using the set_value function']);
            end


            center_ball = obj.center_ball;
            grid_x = state_partition.get_values.partition.grid_x;

            beta = obj.chol_factorization\(obj.chol_factorization'\obj.objective_cost);

            temp_kernel = obj.kernel_function(grid_x',grid_x',obj.kernel_parameter,beta,'Normal');
            
            obj.results_optimisation.opt_obj = max(0,center_ball'*obj.objective_cost ...
                                                  - (obj.radius_ball)*temp_kernel);
        end

        function obj = solve_optimisation_QP(obj,state_partition)
            
            % WE MAY WANT TO DELETE THIS FUNCTION AT A LATER VERSION OF THE
            % CODE. THE NUMERICAL RESULTS WITH THIS DOES NOT LOOK
            % PROMISSING. 

            if isempty(obj.kernel_parameter) || isempty(obj.number_of_samples_KME) ...
                    || isempty(obj.current_state) || isempty(obj.current_input)
                error(['You must first initialize the kernel_parameter and number_of_data fiels of ' ...
                    'the object, as well as current_state and current_input using ' ...
                    'the set_value function']);
            end

            grid_x = state_partition.get_values.partition.grid_x;
            grid_length = length(grid_x); %number of sample points 

            center_ball = sdpvar(grid_length,1);
            
            % Defining constraints
            constraints = [sum(center_ball) == 1, center_ball >= 0]; % alpha is a probability vector
            tempKernel = obj.chol_factorization; tempKernel = tempKernel'*tempKernel;  % Kernel matrix
            constraints = [constraints,center_ball'*tempKernel*center_ball ...
                            - 2/grid_length*center_ball'*tempKernel*ones(grid_length,1) ...
                                + (1/grid_length^2)*ones(grid_length,1)'*tempKernel*ones(grid_length,1) ...
                                            <= obj.radius_ball^2];
            
            options = sdpsettings('solver','mosek','verbose',0); % options to solve the optimisation problem
            lambda = 10;
            objective = obj.objective_cost'*center_ball + lambda*norm(center_ball,1); % objective function

            sol = optimize(constraints,objective,options);

            % Error if unfeasible
            if sol.problem == 1
                error('Unfeasible optimization problem')
            end
            
            % Saving the results of the solution to the optimization
            % problems

            obj.results_optimisation = [];

            obj.results_optimisation.solver_status = sol; % information about opt problem
            obj.results_optimisation.number_variables = length(depends(constraints)); % number of optimisation variable (this needs to be doubled-checked)
            obj.results_optimisation.opt_obj = value(objective) - lambda*value(norm(center_ball,1)); % optimal objective
            
            obj.optimal_distribution = value(center_ball); % saving the result in the public property

        end
        
        function obj = solve_optimisation_KME(obj,state_partition,data_KME)

            obj.results_optimisation = [];

            mod_objective_cost = modify_obj_func(obj.objective_cost,state_partition,...
                data_KME);

            temp = (obj.chol_factorization)\(obj.chol_factorization'...
                *obj.kernel_function(obj.current_state,obj.current_input));

%             obj.kernel_function(obj.current_state,obj.current_input)
%             obj.current_state
%             obj.current_input
%             mod_objective_cost'*temp


            %temp = temp/(abs(temp)); % I AM NOT SURE ABOUT THIS
            temp = abs(temp)/sum(abs(temp)); % I AM NOT SURE ABOUT THIS
            
            %mod_objective_cost'*temp

            %norm_RKHS_mod_objective = sqrt(mod_objective_cost'*obj.chol_factorization'...
            %*obj.chol_factorization*mod_objective_cost)...
            %/norm(obj.chol_factorization'*obj.chol_factorization); % I AM NOT SURE ABOUT THIS

            switch size(data_KME,2)
                case 3
                    obj.results_optimisation.opt_obj = max(0,min(1,obj.eta_param*mod_objective_cost'*temp ...
                                - obj.radius_ball*obj.norm_value_func_RHKS));

                case 5
                    obj.results_optimisation.opt_obj = obj.eta_param*mod_objective_cost'*temp ...
                                                        - obj.radius_ball*obj.norm_value_func_RHKS;
                otherwise
                    not_implemented();
            end
            %1/obj.norm_value_func_RHKS*[mod_objective_cost'*temp,obj.norm_value_func_RHKS]
            %obj.results_optimisation.opt_obj
            %obj.results_optimisation.opt_obj
        end
    end
end


function mod_objective_cost = modify_obj_func(objective_cost,state_partition,data_KME)
    
    number_of_points_KME = size(data_KME,1);
    mod_objective_cost = zeros(number_of_points_KME,1);

    switch size(data_KME,2)
        case 3
            grid_points_KME = data_KME(:,1);

            for i=1:number_of_points_KME
                index_of_point = state_partition.get_element_partition(grid_points_KME(i,:)').index;

                mod_objective_cost(i) = objective_cost(index_of_point);

            end

        case 5
            grid_points_KME = data_KME(:,1:2);
            grid_x = state_partition.partition.grid_x;

            for i=1:number_of_points_KME
                temp = state_partition.get_element_partition(grid_points_KME(i,:)');
                index_of_point = intersect(find(grid_x(:,1) == temp.x(1)),find(grid_x(:,2) == temp.x(2)));
                mod_objective_cost(i) = objective_cost(index_of_point);

            end


        otherwise
            not_implemented();
    end

    
end


