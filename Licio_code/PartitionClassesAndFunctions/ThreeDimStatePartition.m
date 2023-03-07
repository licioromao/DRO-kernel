classdef ThreeDimStatePartition < StatePartition
    
    % This class implements all abstract methods of StatePartition for a
    % three-dimensional grid

    properties
        partition
    end

    methods
        function obj = ThreeDimStatePartition(grid,number_of_points,safe_set, ... 
                                                type_vector_field)
            
            if any(safe_set(:,1) < grid(:,1)) || any(safe_set(:,2) > grid(:,2))
                error('The safe set must be contained inside the grid box');
            end

            if size(grid,1) ~= 3
                error('This class only implements three-dimensional gridding')
            end

            obj = obj@StatePartition(type_vector_field);
            obj.partition = generate_partition(grid,number_of_points);
        end

        function [grid_with_inputs,grid_no_inputs] = create_list(obj,input_partition)

            % This function creates a vector of array that contains the labels of the
            % states of the generated discrete model
            %
            % Inputs: input_partition -- this is the output of the function
            % generatePartition as defined above
            %
            % Ouput: out -- labels for the discrete states
            %
            % This function produces a cell list with the names to each grid element of
            % the partition in the state space. It also creates a cell list for each
            % state-input pair



            number_of_points_axis_1 = size(obj.partition.X1,2);
            number_of_points_axis_2 = size(obj.partition.X2,1);
            number_of_points_axis_3 = size(obj.partition.X3,3);
            number_of_inputs = size(input_partition,1);

            grid_with_inputs = cell(number_of_points_axis_1*number_of_points_axis_2...
                          *number_of_points_axis_3*number_of_inputs,1);
            grid_no_inputs = cell(number_of_points_axis_1*number_of_points_axis_2...
                            *number_of_points_axis_3,1);
            
            % Iterating over the all axis
            for i1 = 1:number_of_points_axis_1
                for i2 = 1:number_of_points_axis_2
                    for i3 = 1:number_of_points_axis_3
                        x = [obj.partition.X1(i2,i1,i3),...
                                obj.partition.X2(i2,i1,i3),...
                                    obj.partition.X3(i2,i1,i3)];

                        % Iterating over the inputs
                        for j=1:number_of_inputs
                            u = input_partition(j,:);

                            temp_index = remaining_iterations(4,[[i1;i2;i3;j],...
                                [number_of_points_axis_1;number_of_points_axis_2;...
                                    number_of_points_axis_3;number_of_inputs]],1,[]);

                            grid_with_inputs{temp_index} = create_x_and_u_string(x,u);
                        end

                        temp_index_2 = remaining_iterations(3,[[i1;i2;i3],...
                            [number_of_points_axis_1;number_of_points_axis_2;...
                                number_of_points_axis_3]],1,[]);

                        grid_no_inputs{temp_index_2} = create_x_and_u_string(x,[]);
                    end
                end
            end
        end
        
        function out = get_size_partition(obj)

            size_axis_1 = obj.partition.X1(1,2,1) - obj.partition.X1(1,1,1);
            size_axis_2 = obj.partition.X2(2,1,1) - obj.partition.X2(1,1,1);
            size_axis_3 = obj.partition.X3(1,1,2) - obj.partition.X3(1,1,1);

            out = [size_axis_1;size_axis_2;size_axis_3];   
        end
        
        function out = get_element_partition(obj,x)
            % This function returns the element of the partition that x belongs to. If
            % the input x is outside the safe set, this function returns an empty
            % array.

            if length(x) ~= 3
                error('The input must be a 3-dimensional vector')
            end

            index = [find(obj.partition.X1(1,:,1) <= x(1),1,'last');
                find(obj.partition.X2(:,1,1) <= x(2),1,'last');
                find(obj.partition.X3(1,1,:) <= x(3),1,'last')]; % Check the element of the partition

            if length(index) ~= 3 % If there is one dimension for which index is empty, return an empty array
                out.index = [1;1;1];
                out.x = obj.partition.grid_x(1,:)';
                out.x = [out.x];
            else
                x_hat = [obj.partition.X1(1,index(1),1);obj.partition.X2(index(2),1,1);obj.partition.X3(1,1,index(3))];
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

                case 'Fishery'

                    VF = VectorFieldFishery(current_state,current_input,param);
                    VF.Noise = noise;
                    next_state = VF.IterateDynamics;

                    temp = obj.get_element_partition(next_state);
                    out.element_partition = temp.index;
                    out.next_state = obj.get_center_partition(temp.x);

                otherwise
                    not_implemented();
            end
        end
   
    end
end