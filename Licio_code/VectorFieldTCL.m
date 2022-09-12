classdef VectorFieldTCL
    
    
    properties
        NextState
        Noise
    end
    
    properties (Access = private)
        CurrentState % stores the current state
        Input % input applied to the dynamics
        R % Thermal resistence in Celsius/kW
        C % Thermal capacity ini kW/Celsius
        P % Range of energy transfer to or from the thermal mass in kW
        eta % Control efficiency
        h % Discretization time in hours
        alpha
        theta
        
        
    end
    
    methods
        function obj = VectorFieldTCL(x,u,param)
            
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.R = param.R;
            obj.C = param.C;
            obj.P = param.P;
            obj.eta = param.eta;
            obj.h = param.h;
            obj.alpha = param.alpha;
            obj.theta = param.theta;
            obj.NextState = [];
            obj.Noise = [];
        end
        
        function obj = SetValues(obj,x,u,param)
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.R = param.R;
            obj.C = param.C;
            obj.P = param.P;
            obj.eta = param.eta;
            obj.h = param.h;
            obj.alpha = param.alpha;
            obj.theta = param.theta;
            
        end
        
        function out = GetValues(obj)
            out.CurrentState = obj.CurrentState;
            out.Input = obj.Input;
            
            out.R = obj.R;
            out.C = obj.C;
            out.P = obj.P;
            out.eta = obj.eta;
            out.h = obj.h;
            out.alpha = obj.alpha;
            out.theta = obj.theta;
            
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
                
                out = obj.CurrentState*obj.alpha + (1-obj.alpha)*(obj.theta - obj.eta*obj.R*obj.P*obj.Input) + obj.Noise;
                
                
            end
            
        end
    end
    
    methods (Access =  private)
    end
end

function VerifyArguments(x,u,param)

if ~isvector(x) || ischar(x)
    error('The first argument must be a row or column vector of numbers');
elseif length(x) ~= 1
    error('The dimension of the first argument must be 1');
end

if ~isvector(u) || ischar(u)
    error('The second argument must be a row or column vector of numbers');
elseif length(u) ~= 1
    error('The dimension of the second argument must be 1');
end

if ~isfield(param,'R') || ~isfield(param,'C') || ~isfield(param,'P') || ~isfield(param,'eta') || ~isfield(param,'h') || ~isfield(param,'alpha') || ~isfield(param,'theta')
    error('The third argument must a structure with fields R, C, P, eta, h, theta, and alpha');
end

end

