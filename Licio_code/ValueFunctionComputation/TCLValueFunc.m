classdef TCLValueFunc < ComputeValueFunction

    properties
        type_vector_field % This value stores the type of Vector field
        param % structure with all the parameters of the problem
        index_safe_set 
    end

    methods
        function obj = TCLValueFunc(number_of_points,time_horizon,type_vector_field,param)
            
            if ~strcmp(type_vector_field,'TCL')
                error('Type of vector field must be TCL')
            end

            obj = obj@ComputeValueFunction(number_of_points,time_horizon);
            
            obj.type_vector_field = type_vector_field;
            obj.param = param;
            obj.index_safe_set = [];
        end

        function obj = get_index_safety(obj,state_partition)
             % Retuns the indices of the safe set associated with the partition
            
            if ~strcmp(state_partition.type_vector_field,'TCL') 
                error('This function can only be used for TCL vector field. Please check if this is your case.')
            end
            
            temp_state_partition = state_partition.partition;
            number_of_points = size(temp_state_partition.grid_x,1);

            % Information about the safe and reach sets
            safe_set = obj.param.safe_set;
            safe_index = false(number_of_points,1);
            
%             switch obj.type_vector_field
%                
%                 case 'ChainInt'
%                     number_of_points = size(temp_state_partition.grid_x,1);
%                 otherwise
%                     NotImplemented();
%             end
 
            % Iterating over states
            for i = 1:number_of_points
                x = temp_state_partition.grid_x(i);  % current state
                
                % if current state belongs to safe set, then set the
                % corresponding index to true
                if min(safe_set(1) <= x) && min(x <= safe_set(2))
                    safe_index(i) = true;
                end
            end  
            % saving the results
            obj.index_safe_set = safe_index;
            
        end

    end
end