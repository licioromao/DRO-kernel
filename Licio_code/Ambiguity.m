classdef (Abstract) Ambiguity
    % This class solves an optimization problem that is relevant to study
    % of distributionally robust ambiguity sets. It is often the case that
    % one needs an optimization problem given by:
    
    %%
    % $\mathrm{minimize} \quad  c^\top p$
    %%
    % $\mathrm{subject~to} \quad p \in \mathcal{A}$
    %%
    % 
    % * $p$ is a n-dimensional vector 
    % * $\mathcal{A}$ is the ambiguity set
    
    % This class contains as public variables p and OptRes, which are the
    % optimal value and minimizer of the optimization problem above. In
    % this implementation, c is a private variable of the class. The
    % definition for all methods that manipulate these variables is given
    % below.
    %%
    
    % Two public variables of the class
    properties
        OptRes 
        p
    end
    
    % Private variable that defines the objective function for the
    % optimization problem above
    properties (Access = protected)
        c         
    end
        
    % Methods that manupulate the variables of the class
    methods
        % Constructor of this class
        function obj = Ambiguity(arg1)
            
            out = obj.verifyArg(arg1);
            obj.c = out;
            obj.OptRes = [];
            p = [];
        end
        
        % Interface with the private parameters of the class
        function obj = setValues(obj,arg1)
            
            out = obj.verifyArg(arg1);
            
            obj.c = out;   % Modifies the private parameter                     
        end
        
        function out = getValues(obj)
            out.c = obj.c; % Gets the value of the private parameter
        end    
    end
    
    % This is a private method that can be invoked by class and subclasses
    methods (Access = protected, Static)
        function out = verifyArg(arg1)
            % This function receives as input a n-dimensional vector arg1
            % and creates a column vector from it. It returns an error is the input
            % is not a n-dimensional vector.
            
            if (isvector(arg1) && ~ischar(arg1)) || isempty(arg1)
                if size(arg1,1) == 1
                    out = arg1';
                else
                    out = arg1;
                end
            else
                error('The private property must be a row or column vector of numbers');
            end  
        end
    end
    
    % Defines the method to solve the optimization problem above as
    % abstract so that we can have different implementation of the method
    methods (Abstract)
        SolveOptimization(obj)
    end
    
end

