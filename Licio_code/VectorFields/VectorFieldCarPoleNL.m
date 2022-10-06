classdef VectorFieldCarPoleNL
    
    
    properties
        NextState
        Noise
    end
    
    properties (Access = private)
        CurrentState % stores the current state
        Input % input applied to the dynamics
        T % Discretization time
        m % mass of the pole
        l % length of the pole
        mt % total mass
        g % gravity
    end
    
    methods
        function obj = VectorFieldCarPoleNL(x,u,param)
            
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.T = param.T;
            obj.m = param.m;
            obj.l = param.l;
            obj.mt = param.mt;
            obj.g = param.g;            
            
            obj.NextState = [];
            obj.Noise = [];
        end
        
        function obj = SetValues(obj,x,u,param)
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.T = param.T;
            obj.m = param.m;
            obj.l = param.l;
            obj.mt = param.mt;
            obj.g = param.g;
            
            
            obj.NextState = [];
            obj.Noise = [];
            
        end
        
        function out = GetValues(obj)
            
            out.CurrentState = obj.CurrentState;
            out.Input = obj.Input;
            
            out.T = obj.T;
            out.m = obj.m;
            out.l = obj.l;
            out.mt = obj.mt;
            out.g = obj.g;
            
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
                tempNonLinearDynamics = @(t,x) NonLinearDynamics(t,x,obj.Input,obj.m,obj.l,obj.g,obj.mt);
                tempOut = ode15s(tempNonLinearDynamics,[0,obj.T],obj.CurrentState);
                out = tempOut.y(:,end);
            end
            
        end
    end
    
end

function VerifyArguments(x,u,param)

if ~isvector(x) || ischar(x)
    error('The first argument must be a row or column vector of numbers');
elseif length(x) ~= 4
    error('The dimension of the first argument must be 1');
end

if ~isvector(u) || ischar(u)
    error('The second argument must be a row or column vector of numbers');
elseif length(u) ~= 1
    error('The dimension of the second argument must be 1');
end

if ~isfield(param,'T') || ~isfield(param,'m') || ~isfield(param,'l') || ...
         ~isfield(param,'mt') || ~isfield(param,'g')
    error('The third argument must a structure with field T, m, l, mt, and g');
end

end

function out = NonLinearDynamics(t,x,u,m,l,g,mt)

% x = [position, velocity, angular position, angular velocity]

out = zeros(4,1);
denominator = l*mt*(4/3 - m*((cos(x(3))^2)/(mt)));

term =  (u+m*l*x(4)^2*sin(x(3)))/(mt);

out(1) = x(2);
out(2) = term  - ...
    (m*l*(g*sin(x(3)) - cos(x(3)))*term*cos(x(3)))/denominator;
out(3) = x(4);
out(4) = (g*sin(x(3)) - cos(x(3))*term*cos(x(3)))/(denominator);

end


