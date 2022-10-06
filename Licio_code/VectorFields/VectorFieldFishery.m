classdef VectorFieldFishery
    
    % This class implements the vector field of the example in Section 6.1
    % of paper "Verification of discrete time stochastic hybrid systems: A
    % stochastic reach-avoid decision problem" by Summers and Lygeros,
    % Automatica, 2010.
    
    properties
        NextState
        Noise
    end
    
    properties (Access = private)
        CurrentState
        Input
        r
        L
        C
    end
    
    methods
        function obj = VectorFieldFishery(x,u,param)
           
            
            VerifyArguments(x,u,param);
                        
            obj.CurrentState = x;
            obj.Input = u;
                  
            obj.r = param.r;
            obj.L = param.L;
            obj.C = param.C;
            obj.NextState = [];
            obj.Noise = [];
        end
        
        function obj = SetValues(obj,x,u,param)
            
            VerifyArguments();
            
            obj.CurrentState = x;
            obj.Input = u;
                 
            obj.r = param.r;
            obj.L = param.L;
            obj.C = param.C;
        end
        
        function out = GetValues(obj)
            out.CurrentState = obj.CurrentState;
            out.Input = obj.Input;
            out.r = obj.r;
            out.L = obj.L;
            out.C = obj.C;
        end
        
        function out = IterateDynamics(obj)
            
            % Equations (14)-(16) of the paper mentioned at the top of this file.
            
            % Input: Noise - a 7-dimensional vector containing the sources of noise in
            %        the model
            %        Input - Control input
            %        CurrentState - current state of the model
            %        param - is a structure that contains all the parameters of the
            %        model
            %
            % Output: out - next state of the model
            
            if isempty(obj.CurrentState) || isempty(obj.Noise) || isempty(obj.Input)
                error('Either private value CurrentState or public property Noise is empty')
            end
                        
            v1 = obj.Noise(1);
            v2 = obj.Noise(2);
            gamma1 = obj.Noise(3);
            gamma2 = obj.Noise(4);
            lambda = obj.Noise(5);
            delta1 = obj.Noise(6);
            delta2 = obj.Noise(7);
            
            R1 = vector_field_R(obj,obj.CurrentState(1));
            R2 = vector_field_R(obj,obj.CurrentState(2));
            
            M = vector_field_M(obj,obj.CurrentState(1),obj.CurrentState(2),obj.Input(1),obj.Input(2),v1,v2,delta1,delta2);
            
            [C1,mode] = vector_field_C(obj,obj.CurrentState(1)+obj.CurrentState(2),obj.Input(1),delta1);
            [C2,~] = vector_field_C(obj,obj.CurrentState(1)+obj.CurrentState(2),obj.Input(2),delta2);
            
            out = [(1-v1)*obj.CurrentState(1) + gamma1*R1 - C1 - lambda*M; % Equation (14)
                (1-v2)*obj.CurrentState(2) + gamma2*R2 - C2 + lambda*M; % Equation (15)
                obj.CurrentState(3) + C1 + C2; % Equation (16)
                mode]; % Mode of the system
        end
    end
    
    methods (Access =  private)
        function out = vector_field_R(obj,x)
            
            % This function implements the function R that composes the vector field,
            % as given in the first column of page 1956 of the paper mentioned at the
            % top of this file.
            
            % Inputs:   x - current state
                       
            out = max(0,obj.r*x*(1-(x/obj.L))); % output of the function
        end
        
        function [out,mode] = vector_field_C(obj,x,d,delta)
            
            % This function implements the function C as shown in page 1956 of the
            % paper at the top of this file.
            
            %  Inputs: x - current state
            %          d - represents the input of the dynamics
            %          delta - is the variability of the catch
            %          param -  is a structure containing L and C as fields
            %
            %  Output: out - value of the expression
            %          mode - mode of the system. Note that this system has two modes,
            %          according to the value of the current state and L
            %
         
            if x < obj.L % if x is less than L
                out = delta*d*obj.C*x/obj.L;
                mode = 1;
            else % otherwise
                out = delta*d*obj.C;
                mode = 2;
            end
        end
        
        function out = vector_field_E(obj,x1,x2,d,v,delta)
            
            % This function implements the function E as shown in page 1956 of the
            % paper at the top of this file.
            
            % Inputs: x1,x2 -  fish biomass in each patch
            %         d - represents the input of the dynamics
            %         delta - is the variability of the catch
            %         v - natural fish mortality
            %         param - is a structure containing L and C as fields
            
            out = (1-v)*x1 - obj.vector_field_C(x1+x2,d,delta); % output of this expression as shown in the paper
            
        end
        
        function out = vector_field_M(obj,x1,x2,d1,d2,v1,v2,delta1,delta2)
            
            % This function implements the function M as shown in page 1956 of the
            % paper at the top of this file, which represents the fish escapament
            % between patch at a given period
            
            % Inputs: x1,x2 -  fish biomass in each patch
            %         d1,d2 - represents the inputs for each patch
            %         delta1,delta2 - is the catch variability in each patch
            %         v1,v2 - natural fish mortality in each patch
            %         param - is a structure containing L and C as fields
            
            
            temp_1 = obj.vector_field_E(x1,x2,d1,v1,delta1);
            temp_2 = obj.vector_field_E(x2,x1,d2,v2,delta2);
            
            out = temp_1 - temp_2;
        end
    end
end

function VerifyArguments(x,u,param)

if ~isvector(x) || ischar(x)
    error('The first argument must be a row or column vector of numbers');
elseif length(x) ~= 3
    error('The dimension of the first argument must be 3');
end

if ~isvector(u) || ischar(u)
    error('The second argument must be a row or column vector of numbers');
elseif length(u) ~= 2
    error('The dimension of the second argument must be 2');
end

if ~isfield(param,'r') || ~isfield(param,'L') || ~isfield(param,'C')
    error('The third argument must a structure with fields r, L, and C');
end

end

