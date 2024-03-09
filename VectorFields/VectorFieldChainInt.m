classdef VectorFieldChainInt
    
    % This vector field has been studied in Section 5.1 in the paper
    
    properties
        NextState
        Noise
    end
    
    properties (Access = private)
        CurrentState
        Input
        T
    end
    
    methods
        
        function obj = VectorFieldChainInt(x,u,param)
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.T = param.T;
            
            obj.NextState = [];
            obj.Noise = [];
            
        end
        
        function obj = SetValues(obj,x,u,param)
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.T = param.T;
            
            obj.NextState = [];
            obj.Noise = [];
            
        end
        
        function out = GetValues(obj)
            
            out.CurrentState = obj.CurrentState;
            out.Input = obj.Input;
            
            out.T = obj.T;
            
        end
        
         function out = IterateDynamics(obj)
            
            % This is the vector field of Insoon's paper called "A dynamic game approach to distributionally robust safety specifications for stochastic systems" in Automatica.
            
            % Input: Noise - This is the noise affecting the system
            %        Input - Control input
            %        CurrentState - current state of the model
            %        param - is a structure that contains all the parameters of the
            %        model
            %
            % Output: out - next state of the model
            if isempty(obj.Noise)
                error('You must first initialize the field Noise');
            else
                
                out = [1,obj.T;0,1]*obj.CurrentState + [obj.T^2/2;obj.T]*obj.Input + obj.Noise;
            end
            
        end
        
    end
    
    
end

function VerifyArguments(x,u,param)

if ~isvector(x) || ischar(x)
    error('The first argument must be a row or column vector of numbers');
elseif length(x) ~= 2
    error('The dimension of the first argument must be 1');
end

if ~isvector(u) || ischar(u)
    error('The second argument must be a row or column vector of numbers');
elseif length(u) ~= 1
    error('The dimension of the second argument must be 1');
end

if ~isfield(param,'T') 
    error('The third argument must a structure with field T');
end

end