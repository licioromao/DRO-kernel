classdef KernelBasedAmbiguity < DistanceBasedAmbiguity
    % This class implements the computation of the Kernel-based ambiguity
    % set problem. It inherits from the class DistanceBasedAmbiguity
    
    % This is very specific of the Gaussian Kernel. Can be generalized later.    
    properties
        param % This variable must be made global at some point
        
        number_of_data
        
        current_state
        current_input
        
        kernel_function % This parameter is a handle to a function that defined the Kernel on the state space
        kernel_parameter
        chol_factorization
    end
       
    methods
        function obj = KernelBasedAmbiguity(objective_cost,radius_ball, ... 
                        center_ball,kernel_function,kernel_parameter,grid_axis)

            obj = obj@DistanceBasedAmbiguity(objective_cost,radius_ball,center_ball); % Calling the connstructor of the parent class
                
            obj.param = []; % This should be made a global variable at some point
            
            obj.number_of_data = [];
           
            obj.current_state = [];
            obj.current_input = [];
            
            % Assigning the kernel parameters to the object
            obj.kernel_parameter = kernel_parameter;
            obj.kernel_function = kernel_function;
            obj.chol_factorization = kernel_function(grid_axis',grid_axis',kernel_parameter,[],'MatrixFac');
        end
                       
        function obj = solve_optimisation_conservative(obj,obj_Partition)
            % This function implements the solution of a kernel ambiguity
            % minimisation where the center of the ambiguity ball is equal
            % to the empirical estimate. 
            % WE MAY WANT TO DELETE THIS FUNCTION AT A LATER VERSION OF THE
            % CODE. THE NUMERICAL RESULTS WITH THIS DOES NOT LOOK
            % PROMISSING. 

            if isempty(obj.kernel_parameter) || isempty(obj.number_of_data) ... 
                      || isempty(obj.current_state) || isempty(obj.current_input)
                error(['You must first initialize the kernel_parameter and number_of_data ' ...
                    'fiels of the object, as well as current_state and current_input ' ...
                    'using the set_value method']);
            end

            center_ball = obj.center_ball;
            grid_axis = obj_Partition.get_values.partition.grid_axis; % CHECK THIS AFTER CHANGING PARTITION OBJECT

            tempKernel = obj.kernel_function(grid_axis',grid_axis', ... 
                                    obj.kernel_parameter,center_ball,'Normal');
            

            obj.results_optimisation.opt_obj = max(0,center_ball'*obj.objective_cost -... 
                                        (obj.radius_ball)*tempKernel-2*obj.radius_ball);
        end

        function obj = solve_optimisation(obj,obj_Partition)

            % WE MAY WANT TO DELETE THIS FUNCTION AT A LATER VERSION OF THE
            % CODE. THE NUMERICAL RESULTS WITH THIS DOES NOT LOOK
            % PROMISSING. 
            if isempty(obj.kernel_parameter) || isempty(obj.number_of_data) ...
                  || isempty(obj.current_state) || isempty(obj.current_input)
                error(['You must first initialize the kernel_parameter and ' ...
                    'number_of_data fiels of the object, as well as ' ...
                    'current_state and current_input using the set_value function']);
            end


            center_ball = obj.center_ball;
            grid_axis = obj_Partition.get_values.partition.grid_axis;

            beta = obj.chol_factorization\(obj.chol_factorization'\obj.objective_cost);

            temp_kernel = obj.kernel_function(grid_axis',grid_axis',obj.kernel_parameter,beta,'Normal');
            
            obj.results_optimisation.opt_obj = max(0,center_ball'*obj.objective_cost ...
                                                  - (obj.radius_ball)*temp_kernel);
        end

        function obj = solve_optimisation_QP(obj,obj_Partition)
            
            % WE MAY WANT TO DELETE THIS FUNCTION AT A LATER VERSION OF THE
            % CODE. THE NUMERICAL RESULTS WITH THIS DOES NOT LOOK
            % PROMISSING. 

            if isempty(obj.kernel_parameter) || isempty(obj.number_of_data) ...
                    || isempty(obj.current_state) || isempty(obj.current_input)
                error(['You must first initialize the kernel_parameter and number_of_data fiels of ' ...
                    'the object, as well as current_state and current_input using ' ...
                    'the set_value function']);
            end

            grid_axis = obj_Partition.get_values.partition.grid_axis;
            grid_length = length(grid_axis); %number of sample points 

            center_ball = sdpvar(grid_length,1);
            
            % Defining constraints
            constraints = [sum(center_ball) == 1, center_ball >= 0]; % alpha is a probability vector
            tempKernel = obj.chol_factorization; tempKernel = tempKernel'*tempKernel;  % Kernel matrix
            constraints = [constraints,center_ball'*tempKernel*center_ball ...
                            - 2/grid_length*center_ball'*tempKernel*ones(grid_length,1) ...
                                + (1/grid_length^2)*ones(grid_length,1)'*tempKernel*ones(grid_length,1) ...
                                            <= obj.epsilon^2];
            
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
            obj.results_optimisation.optimal_obj = value(objective) - lambda*value(norm(center_ball,1)); % optimal objective
            
            obj.optimal_distribution = value(center_ball); % saving the result in the public property

        end
        
    end
end


