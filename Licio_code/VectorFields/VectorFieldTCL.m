classdef VectorFieldTCL
    
    
    properties
        next_state
        noise
    end
    
    properties (Access = private)
        current_state % stores the current state
        input % input applied to the dynamics
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
            
            
            verify_arguments(x,u,param);
            
            obj.current_state = x;
            obj.input = u;
            
            obj.R = param.R;
            obj.C = param.C;
            obj.P = param.P;
            obj.eta = param.eta;
            obj.h = param.h;
            obj.alpha = param.alpha;
            obj.theta = param.theta;
            obj.next_state = [];
            obj.noise = [];
        end
        
        function obj = set_values(obj,x,u,param)
            
            verify_arguments(x,u,param);
            
            obj.current_state = x;
            obj.input = u;
            
            obj.R = param.R;
            obj.C = param.C;
            obj.P = param.P;
            obj.eta = param.eta;
            obj.h = param.h;
            obj.alpha = param.alpha;
            obj.theta = param.theta;
            
        end
        
        function out = get_values(obj)

            out.current_state = obj.current_state;
            out.input = obj.input;
            
            out.R = obj.R;
            out.C = obj.C;
            out.P = obj.P;
            out.eta = obj.eta;
            out.h = obj.h;
            out.alpha = obj.alpha;
            out.theta = obj.theta;
            
        end
        
        function out = iterate_dynamics(obj)
            
            % This is the vector field of Insoon's paper called "A dynamic game approach to distributionally robust safety specifications for stochastic systems" in Automatica.
            
            % Input: noise - This is the noise affecting the system
            %        input - Control input
            %        current_state - current state of the model
            %        param - is a structure that contains all the parameters of the
            %        model
            %
            % Output: out - next state of the model
            if isempty(obj.noise)
                error('You must first initialize the field Noise');
            else
                
                out = obj.current_state*obj.alpha + (1-obj.alpha)*(obj.theta...
                                      -obj.eta*obj.R*obj.P*obj.input) + obj.noise;  
            end
            
        end
    end
    
end

function verify_arguments(x,u,param)

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

if ~isfield(param,'R') || ~isfield(param,'C') || ~isfield(param,'P') ||...
                ~isfield(param,'eta') || ~isfield(param,'h') || ~isfield(param,'alpha')...
                                                || ~isfield(param,'theta')
    error('The third argument must a structure with fields R, C, P, eta, h, theta, and alpha');
end

end

