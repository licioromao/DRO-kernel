classdef OneDimStatePartition < StatePartition
    
    % This class implements all abstract methods of StatePartition for a
    % one-dimensional grid

    properties
        partition
    end

    methods
        function obj = OneDimStatePartition(grid,number_of_points,safe_set, ... 
                                                type_vector_field)
            
            if safe_set(1) < grid(1) || safe_set(2) > grid(2)
                error('The safe set must be contained inside the grid box');
            end

            if size(grid,1) ~= 1
                error('This class only implements one-dimensional gridding');
            end

            obj = obj@StatePartition(type_vector_field);
            obj.partition = generate_partition(grid,number_of_points);
        end

        function [grid_with_inputs,grid_no_inputs] = create_list(obj,input_partition)

            %             This function creates a vector of array that contains the labels of the
            %             states of the generated discrete model
            %
            %             Inputs: input_partition -- this is the output of the function
            %             generatePartition as defined above
            %
            %             Ouput: out -- labels for the discrete states
            %
            % This function produces a cell list with the names to each grid element of
            % the partition in the state space. It also creates a cell list for each
            % state-input pair
            
            state_partition = obj.partition;

            number_of_points = size(state_partition.X,1);
            number_of_inputs = size(input_partition,1);

            grid_with_inputs = cell(number_of_inputs*number_of_points,1);
            grid_no_inputs = cell(number_of_points,1);

            % Iterating over the grid points
            for i = 1:number_of_points
                x = state_partition.X(i);

                % Iterating over the inputs
                for j=1:number_of_inputs

                    temp_index = remaining_iterations(2,[[i;j],[number_of_points;number_of_inputs]],1,[]);
                    grid_with_inputs{temp_index} = create_x_and_u_string(x,input_partition(j));
                end
                grid_no_inputs{i} = create_x_and_u_string(x,[]);
            end
        end
        
        function out = get_values(obj)
            out.type_vector_field = obj.type_vector_field;
            out.partition = obj.partition;
        end

        function out = get_size_partition(obj)
            out = obj.partition.X(2) - obj.partition.X(1);
        end
        
        function out = get_element_partition(obj,x)
            % This function returns the element of the partition that x belongs to. If
            % the input x is outside the safe set, this function returns an empty
            % array.

            index = find(obj.partition.X <= x,1,'last');

            if isempty(index)
                out.index = 1;
                out.x = obj.partition.X(1);
            else
                x_hat = obj.partition.X(index);
                out.index = index;
                out.x = x_hat;
            end
        end
        
        function out = get_center_partition(obj,x)
            % This function receives as input a point x in the lattice and returns the
            % center of the partition.

            out = x + obj.get_size_partition()/2;
        end
        
        function out = compute_element_partition(obj,current_state,... 
                current_input,noise,param)

            % Uses the vector field defined by the type_of_vector_field
            % parameter and iterate the dynamics

            out.next_state = [];
            out.element_partition = [];

            switch obj.type_vector_field % This field is inherited from the parent class

                case 'TCL'

                    VF = VectorFieldTCL(current_state,current_input,param);
                    VF.noise = noise;
                    next_state = VF.iterate_dynamics;

                    temp = obj.get_element_partition(next_state);
                    out.element_partition = temp.index;
                    out.next_state = obj.get_center_partition(temp.x);

                otherwise
                    not_implemented();
            end
        end
   
    end
end
