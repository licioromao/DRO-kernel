classdef VectorFieldLTI
    
    
    properties
        next_state
        noise
        current_state % stores the current state
        input % input applied to the dynamics
        A
        B
    end
    
    methods
        function obj = VectorFieldLTI(x,u,A,B)
            
            
            verify_arguments(x,A,B);
            
            obj.current_state = x;
            obj.input = u;

            obj.A = A;
            obj.B = B;
            
            
            obj.next_state = [];
            obj.noise = [];
        end
        
        function obj = set_values(obj,x,u,A,B)
            
            verify_arguments(x,A,B);
            
            obj.current_state = x;
            obj.input = u;
            
            obj.A = A;
            obj.B = B;
            
        end
        
        function out = get_values(obj)

            out.current_state = obj.current_state;
            out.input = obj.input;
            
            out.A = obj.A;
            out.B = obj.B;
            
        end
        
        function out = iterate_dynamics(obj)
            
            % This is a LTI vector field
            
            % Input: noise - This is the noise affecting the system
            %        input - Control input
            %        current_state - current state of the model
            %
            % Output: out - next state of the model

            if isempty(obj.noise)
                error('You must first initialize the field Noise');
            else
                if isempty(obj.B)
                    out = obj.A*obj.current_state + obj.noise;  
                else
                    if size(obj.input,1) ~= size(obj.B,2)
                        error('Number of inputs must be equal to the number of columns of B');
                    end
                    out = obj.A*obj.current_state + obj.B*obj.input + obj.noise;
                end
            end
            
        end
    end
    
end

function [] = verify_arguments(x,A,B)
    
if ~isempty(B)
    if size(A,1) ~= size(A,2)
        error('Matrix A must be square');
    end

    if size(A,1) ~= size(B,1)
        error('The number of rows of A and B must coincide');
    end

    if size(x,1) ~= length(A)
        error('Dimension of matrix A must be equal to the dimension of x');
    end

end

end

