% Function that computes the value function

function value_func = main_value_function_iteration(state_partition,...
                                input_partition,type_vector_field,...
                                    struct_ambiguity_types,bool_signal,param)

if ~bool_signal
    switch struct_ambiguity_types.name
        case 'NoAmbiguity'
            fprintf('\nComputing value function (without ambiguity)...\n')
        case 'MomentAmbiguity'
            fprintf('\nComputing value function (moment ambiguity)...\n')

            switch type_vector_field
                case 'Fishery'
                    value_func = value_func.getIndexReachAvoid(state_partition.getValues);
                case 'TCL'
                    number_of_points = param.number_of_points;
                    time_horizon = param.time_horizon;
                    radius_mean = struct_ambiguity_types.radius_mean;
                    radius_variance = struct_ambiguity_types.radius_variance;


                    value_func = MomentTCLValueFunc(number_of_points,...
                                        time_horizon,type_vector_field,...
                                            radius_mean,radius_variance,...
                                                    param);

                    value_func = value_func.get_index_safety(state_partition.get_values);
                    value_func = value_func.backward_iteration(state_partition,...
                                                                input_partition,...
                                                                param.outer_loop_info);
                case 'ChainInt'
                    value_func = value_func.getIndexSafety(state_partition.getValues);
                case 'CarPole'
                    value_func = value_func.getIndexReachAvoid(state_partition.getValues);
                case 'CarPoleNL'
                    value_func = value_func.getIndexReachAvoid(state_partition.getValues);
                otherwise
                    NotImplemented();
            end


        case 'KernelAmbiguity'
            fprintf('\nComputing value function (kernel ambiguity)...\n')
        case 'KLdivAmbiguity'
            fprintf('\nComputing value function (KL ambiguity)...\n')
        case 'WassersteinAmbiguity'
            fprintf('\nComputing value function (Wasserstein ambiguity)...\n')
        otherwise
            error('Type of Ambiguity not implemented')
    end
    
    
    fprintf('Done\n');
else
    value_func = [];
    switch struct_ambiguity_types.name
        case 'NoAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncNoAmbiguity');
        case 'MomentAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncMoment');
        case 'KernelAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncKernel');
        case 'KLdivAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncKLdivAmbiguity');
        case 'WassersteinAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncWasserstein');
        otherwise
            error('Type of Ambiguity not implemented')
    end
end

end