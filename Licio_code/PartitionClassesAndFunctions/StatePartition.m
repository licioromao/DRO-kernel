classdef StatePartition
       
    properties 
        type_vector_field
    end
    
    methods
        function obj = StatePartition(type_vector_field)
            obj.type_vector_field = type_vector_field;
        end
    end

    % Abstract methods to be implemented by its children
    methods (Abstract)
        get_values(obj)
        create_list(obj,input_partition)
        get_size_partition(obj)
        get_element_partition(obj,x)
        get_center_partition(obj,x)
        compute_element_partition(obj,current_state,current_input,...
            noise,param)
    end
end





