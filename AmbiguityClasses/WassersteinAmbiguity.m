classdef WassersteinAmbiguity < DistanceBasedAmbiguity
    % This subclass solves the distributionally robust optimization problem
    % defined in the parent class Ambiguity (which the parent of
    % DistanceBasedAmbiguity) using the 1-Wasserstein distance. This
    % class, in addition to the properties inherited from its parents,
    % implements the abstract function solve_optimisation, and its thus the
    % first class in the hierarchy that is not abstract
    
    % The other method in the class computes the distance between the
    % probability measure that defines the center of the ambiguity set and
    % in the public field given in p
        
    methods
        function obj = WassersteinAmbiguity(arg1,arg2,arg3)
            obj = obj@DistanceBasedAmbiguity(arg1,arg2,arg3); % Calling the constructor of its parent
        end
        
        function obj = solve_optimisation(obj)
            
            % This function implements the solution of the optimization
            % problem defined in the parent using the 1-Wasserstein metric
            % as defined in the paper "Robust Markov decision Processes
            % with data-driven, distance-based ambiguity sets", pag. 1007,
            % by Sivaramakrishnan Ramani, et all, 2022, SIAM Journal of
            % Optimization. Our implementation is given by equation (4.7)
            % in the aforementioned paper.
            
            m = length(obj.objective_cost); % dimension of the optimization problem
            x = sdpvar(m,m,'full'); % this represents a matrix whose (i,j)-entry is the variable x_{ij} in the paper
            p = sdpvar(m,1); % this is the optimization varible, q in the formulation of the paper is the center of the ambiguity set
            
            tempMatrix = abs(kron(ones(m,1),1:1:m) - kron(ones(m,1),1:1:m)');
            tempMatrix = tempMatrix.*x; % the sum of rows and columns gives the expression sum_{over i and j} of |i-j|*x_{i,j}
            
            % defining the contraints as in the formulation of the paper
            constraints = [];
            
            % These constraints corresponds to an ambiguity set of radius
            % epsilon in the Wasserstein ball
            constraints = [constraints,sum(x,2) == obj.center_ball, sum(x,1)' == p];
            constraints = [constraints,sum(sum(tempMatrix)) <= obj.radius_ball,x>=0];
            constraints = [constraints,sum(p) == 1,p>= 0];
            
            % Setting solver configuration. One needs to have mosek
            % installed
            options = sdpsettings('solver','mosek','verbose',0);
            
            % Solving the problem
            cost_opt = obj.objective_cost'*p; % This is the objective function in the main optimization problem
            sol = optimize(constraints,cost_opt,options);
            
            if sol.problem == 1
                error('Unfeasible optimization problem') % error if unfeasible
            end
           
            % Saving the results of the solution to the optimization
            % problems

            obj.results_optimisation = [];

            obj.results_optimisation.solver_status = sol; % information about opt problem
            obj.results_optimisation.number_variables = length(depends(constraints)); % number of optimisation variable (this needs to be doubled-checked)
            obj.results_optimisation.optimal_obj = value(sum(sum(tempMatrix))); % optimal distance
            obj.results_optimisation.optimal_coupling = value(x); % optimal joint distribution
            
            obj.optimal_distribution = value(p); % saving the result in the public property
        end
        
        function obj = compute_wasserstein_distance(obj)
            
            % This function computes the Wasserstein distance between the
            % center of the ambiguity set (private property q) and the
            % probability distribution given in the public field p
            
            m = length(obj.c); % dimension of the optimization problem
            x = sdpvar(m,m,'full'); % this represents a matrix whose (i,j)-entry is the variable x_{ij} in the paper
            
            % If there is no public field p, then return an error
            if isempty(obj.optimal_distribution)
                error('Public field optimal_distribution is empty. Either assign a valid probability vector or run the function solve_optimisation');
            end
                
            optimal_distribution = obj.optimal_distribution; % creating a local variable for p
            
            tempMatrix = abs(kron(ones(m,1),1:1:m) - kron(ones(m,1),1:1:m)');
            tempMatrix = tempMatrix.*x; % the sum of rows and columns gives the expression sum_{over i and j} of |i-j|*x_{i,j}
            
            constraints = [];
            
            constraints = [constraints,sum(x,2) == obj.center_ball, sum(x,1)' == optimal_distribution, x>=0]; % constraint defining the joint distribution
                        
            options = sdpsettings('solver','mosek','verbose',0);
            
            objective_function = sum(sum(tempMatrix)); % this is the 1-Wasserstein distance between p and q
            sol = optimize(constraints,objective_function,options); % solving the problem
            
            % Error if unfeasible
            if sol.problem == 1
                error('Unfeasible optimization problem')
            end
          
            obj.results_optimisation = [];
            obj.results_optimisation.solver_status = sol; % information about opt problem
            obj.results_optimisation.number_variables = length(depends(constraints)); % number of optimisation variable (this needs to be doubled-checked)
            obj.results_optimisation.optimal_obj = value(sum(sum(tempMatrix))); % distance between the distributions 
            obj.results_optimisation.optimal_coupling = value(x); % optimal joint distrirbution (or optimal coupling)
        end
    end
end

