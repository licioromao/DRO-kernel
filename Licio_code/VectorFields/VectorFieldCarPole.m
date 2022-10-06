classdef VectorFieldCarPole
    
    
    properties
        NextState
        Noise
    end
    
    properties (Access = private)
        CurrentState % stores the current state
        Input % input applied to the dynamics
        T % Discretization time
        A % Continuous-time dynamic matrix
        B % Continuous-time input matrix
        Ad % Discrete-time dynamic matrix
        Bd % Discrete-time input matrix
    end
    
    methods
        function obj = VectorFieldCarPole(x,u,param)
            
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.T = param.T;
            obj.A = param.A;
            obj.B = param.B;
            
            [obj.Ad,obj.Bd] = DiscretizeSystem(obj.A,obj.B,obj.T);
            
            obj.NextState = [];
            obj.Noise = [];
        end
        
        function obj = SetValues(obj,x,u,param)
            
            VerifyArguments(x,u,param);
            
            obj.CurrentState = x;
            obj.Input = u;
            
            obj.T = param.T;
            obj.A = param.A;
            obj.B = param.B;
            
            [obj.Ad,obj.Bd] = DiscretizeSystem(obj.A,obj.B,obj.T);
            
            obj.NextState = [];
            obj.Noise = [];
            
        end
        
        function out = GetValues(obj)
            out.CurrentState = obj.CurrentState;
            out.Input = obj.Input;
            
            out.T = obj.T;
            out.A = obj.A;
            out.B = obj.B;
            out.Ad = obj.Ad;
            out.Bd = obj.Bd;
            
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
                out = obj.Ad*obj.CurrentState + obj.Bd*obj.Input + obj.Noise;
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

if ~isfield(param,'T') || ~isfield(param,'A') || ~isfield(param,'B')
    error('The third argument must a structure with field T, A, B');
end

end

function [Ad,Bd] = DiscretizeSystem(A,B,T)
    
n = size(A,1);

if n ~= size(A,2)
    error('The dynamical matrix must be square');
end

if n~= size(B,1)
    error('Matrix A and B must have the same number of rows');
end

Ad = expm(A*T);

TempFunc = @(u) expm(A*(T-u))*B;

Bd = integral(TempFunc,0,T,'ArrayValued',true);

end

