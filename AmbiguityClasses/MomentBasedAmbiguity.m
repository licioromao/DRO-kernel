classdef MomentBasedAmbiguity < Ambiguity
    % This class is a subclass of the abstract class ambiguity and
    % implements the abstract method SolveOptimization using a moment-based
    % formulation as in the Automatica paper by Insoon Yang "A dynamic game
    % approach to distributionally robust safety specifications for
    % stochastic systems". To define the ambiguity set using moment
    % information we follow the formulation of equation (5) in Insoon's
    % paper with slightly different notation.
    
    
    %  Comparison of the fields of the class with the notation of Insoon's
    %  paper.
    %
    %   * radius_variance is the parameter defining the variance of the ambiguity
    %   set. It stands for the parameter c in equation (5) of Insoon's
    %   paper.
    %
    %   * radius_mean is the parameter defining the expectation radius in Insoon's
    %   paper. It stands for the parameter b in equation (5) of Insoon's
    %   paper.
    %
    %
       
    % Private fields of the class
    properties (Access = protected)
        mean_center  % "center" of the mean
        variance_center % "center" of the variance
        radius_variance % Radius of the "variance" characterizing the ambiguity set
        radius_mean % Radius of the "expectation" characterizing the ambiguity set
        support_set_distribution % This is the discrete set from where the random variable takes value
    end
    
    methods
        % Constructor of the class
        function obj = MomentBasedAmbiguity(objective_cost,covariance_matrix, ...
                         mean_center,radius_variance,radius_mean, ...
                            support_set_distribution)
          
            check_PSD_and_symmetric(covariance_matrix);
            
            obj = obj@Ambiguity(objective_cost); % Calls the constructor of the class ambiguity
            
            % Saving the values of the private field of the class
            obj.support_set_distribution = support_set_distribution;

            obj.mean_center = mean_center;
            obj.variance_center = covariance_matrix;
            obj.radius_mean = radius_mean;
            obj.radius_variance = radius_variance;       
            
        end
        
        function out = get_values(obj)
            
            % This function gets the values of the private fields and
            % returns them as a structure
            
            out.objective_cost = obj.objective_cost;
            out.mean_center = obj.mean_center;
            out.variance_center = obj.variance_center;
            out.radius_variance = obj.radius_variance;
            out.radius_mean = obj.radius_mean;
            out.support_set_distribution = obj.support_set_distribution;
        end
        
        function obj = set_values(obj,objective_cost,covariance_matrix, ...
                                    mean_center,radius_variance,radius_mean,... 
                                        support_set_distribution)
            
            % This function modifies the private properties of the class.
            % It mimics the structure of its constructor 
            
            check_PSD_and_symmetric(covariance_matrix);
                        
            % Updating the values of the private fields
            obj.objective_cost = objective_cost;
            
            obj.support_set_distribution = support_set_distribution;
            obj.mean_center = mean_center;
            obj.variance_center = covariance_matrix;
            
            obj.radius_mean = radius_mean;
            obj.radius_variance = radius_variance; 
            
        end
        
        function obj = solve_optimisation(obj)
            
            % This function solves the optimization problem defined in the
            % description of the class Ambiguity. The formulation below is
            % based on the paper by Insoon Yang mentioned at the top of
            % this file
            
            m = size(obj.support_set_distribution,2); % number of samples
            n = size(obj.support_set_distribution,1);
            
            x = sdpvar(m,1); % optimization variable -> this will result in the optimal distribution
                
            constraints = [];
            constraints = [constraints, x>=0, sum(x) == 1]; % setting the constraint for x to be a probability distribution
            % The next two constraints define the ambiguity set based on
            % moments
            constraints = [constraints, abs(obj.support_set_distribution*x - obj.mean_center)  <= obj.radius_mean];
            temp = (obj.support_set_distribution - obj.mean_center);
            constraints = [constraints, temp*diag(x)*temp'<= obj.radius_variance*obj.variance_center*eye(n)];
            
            options = sdpsettings('solver','mosek','verbose',0); % setting properties to solve the optimization variable
            
            sol = optimize(constraints,obj.objective_cost'*x,options);
            
            % Error if unfeasible
            if sol.problem == 1
                error('Unfeasible optimization problem')
            end
            

            % Saving the results of the solution to the optimization
            % problems

            obj.results_optimisation = [];

            obj.results_optimisation.solver_status = sol; % information about opt problem
            obj.results_optimisation.number_variables = length(depends(constraints)); % number of optimisation variable (this needs to be doubled-checked)
            obj.results_optimisation.optimal_obj = value(obj.objective_cost'*x); % optimal distance
            
            obj.optimal_distribution = value(x); % saving the result in the public property

        end
   
        function obj = solve_optimisation_dual(obj)
            
            % This function is similar to the methods solve_optimisation,
            % with the difference being that it finds the optimal solution
            % by solving the dual problem as described in Theorem 2 in
            % Insoon's paper.
            
            m = size(obj.support_set_distribution,2); % Number of samples
            n = size(obj.support_set_distribution,1); % dimension of the Euclidean space from which the random variable takes value
            
            % These are optimization variable according to the formulation
            % in the paper of Insoon
            lambda_bar = sdpvar(n,1); 
            bar_lambda = sdpvar(n,1);
            Lambda = sdpvar(n,n);
            nu = sdpvar(1);
                
            constraints = [];
            
            constraints = [constraints,bar_lambda >= 0, lambda_bar >=0, Lambda >= 0];
            
            for i = 1:m
                constraints = [constraints, obj.c(i) + (obj.support_set_distribution(:,i) - obj.mean_center)'*Lambda*(obj.support_set_distribution(:,i) - obj.mean_center) + (obj.support_set_distribution(:,i) - obj.mean_center)'*(lambda_bar - bar_lambda) + nu >= 0];
            end
            
            options = sdpsettings('solver','mosek','verbose',0);
            
            objective = bar_lambda'*ones(n,1)*obj.radius_mean + lambda_bar'*ones(n,1)*obj.radius_mean + nu + obj.radius_variance*trace(Lambda*obj.variance_center);
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
            obj.results_optimisation.optimal_obj = value(-objective); % optimal distance
            obj.optimal_distribution = []; % In this formulation, we cannot recover the optimal distribution
        end
    end
end

function out = check_PSD_and_symmetric(covariance_matrix)

tol = 1e-5;

if norm(covariance_matrix - covariance_matrix') < tol && min(real(eig(covariance_matrix))) > -tol ...
            && (ismatrix(covariance_matrix) && (size(covariance_matrix,1) == size(covariance_matrix,2)))
    out = true;
else
    error('The matrix must be square, symmetric and PSD')
end

end




