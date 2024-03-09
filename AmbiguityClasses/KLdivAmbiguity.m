classdef KLdivAmbiguity < DistanceBasedAmbiguity
    % This subclass solves the distributionally robust optimization problem
    % defined in the parent class Ambiguity (which the parent of
    % DistanceBasedAmbiguity) using the KL divergence. This
    % function, in addition to the properties inherited from its parents,
    % implements the abstract function SolveOptimization, and its thus the
    % first class in the hierarchy that is not abstract
    
    % The other method in the class computes the KL divergence between the
    % probability measure that defines the center of the ambiguity set and
    % in the public field given in p
    
    methods
        function obj = KLdivAmbiguity(arg1,arg2,arg3)
            obj = obj@DistanceBasedAmbiguity(arg1,arg2,arg3) % Calling the constructor of its parent
        end
        
        
        function obj = SolveOptimization(obj)
            
            % This function implements the solution of the optimization
            % problem defined in the parent using the KL metric in the
            % space of probability distributions. To compute the KL
            % distance we use the built-in MATLAB function kullbackleibler.
            
            m = length(obj.c); % dimension of the optimization problem
            p = sdpvar(m,1); % this is the optimization varible, q in the formulation of the paper is the center of the ambiguity set
            
            Constraints = [];
            
            % Constraints defining the ambiguity set
            Constraints = [Constraints,kullbackleibler(p,obj.q) <= obj.epsilon];
            Constraints = [Constraints,sum(p) == 1,p>= 0];
            
            options = sdpsettings('solver','mosek','verbose',0); % 
            
            
            % Setting solver configuration. One needs to have mosek
            % installed
            ObjFunc = obj.c'*p; % This is the objective function in the main optimization problem
            sol = optimize(Constraints,ObjFunc,options);
            
            if sol.problem == 1
                error('Unfeasible optimization problem') % error if unfeasible
            end
            
            % Saving the results
            tol = 1e-2;
            
            % All the lines below modifies the public parameter OptRes
            % using the information obtained from the solver
            obj.OptRes = [];
            obj.OptRes.SolverStatus = sol; % information about opt problem
            %obj.OptRes.dualVar = dual(Constraints);
            
            obj.OptRes.Nvar = length(depends(Constraints)); % number of optimisation variable (this needs to be doubled-checked)
            
            obj.OptRes.p = value(p);  % Optimal distribution
            obj.OptRes.opt_value = value(ObjFunc); % optimal value
            
            obj.p = value(p); % Optimal distribution
        end
        
        function obj = KLDistance(obj)
            
            % This function computes the KL divergence between the
            % center of the ambiguity set (private property q) and the
            % probability distribution given in the public field p
            % If there is no public field p, then return an error
            
            if isempty(obj.p)
                error('Public field p is empty. Either assign a valid probability vector or run the function SolveOptimization');
            end
                          
            obj.OptRes = [];
            
            obj.OptRes.opt_value = kullbackleibler(obj.p,obj.q); % Computing the KL divergence distance  
        end
        
    end
end

