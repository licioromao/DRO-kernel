% Function that computes the value function

function value_func = main_value_function_iteration(state_partition,...
                                input_partition,type_vector_field,...
                                    struct_ambiguity_types,bool_signal,param)

if ~bool_signal
    switch type_vector_field
        case 'TCL'
            switch struct_ambiguity_types.name
                case 'NoAmbiguity'
                    fprintf('\nComputing value function (without ambiguity)...\n');

                    number_of_points = param.number_of_points;
                    time_horizon = param.time_horizon;

                    value_func = NoAmbiguityTCLValueFunc(number_of_points,...
                                        time_horizon,type_vector_field,param);

                case 'MomentAmbiguity'
                    fprintf('\nComputing value function (moment ambiguity)...\n')

                    number_of_points = param.number_of_points;
                    time_horizon = param.time_horizon;
                    radius_mean = struct_ambiguity_types.radius_mean;
                    radius_variance = struct_ambiguity_types.radius_variance;


                    value_func = MomentTCLValueFunc(number_of_points,...
                                        time_horizon,type_vector_field,...
                                            radius_mean,radius_variance,...
                                                    param);
                case 'KernelAmbiguity'
                    fprintf('\nComputing value function (kernel ambiguity)...\n');

                    number_of_points = param.number_of_points;
                    time_horizon = param.time_horizon;

                    radius_ball = struct_ambiguity_types.radius_ball;
                    kernel_parameter = struct_ambiguity_types.kernel_parameter;
                    type_value_func_computation = struct_ambiguity_types.type_value_func_computation;

                    if strcmp(type_value_func_computation,'KME')
                        value_func = KernelTCLValueFunc(number_of_points,time_horizon,...
                                                          type_vector_field,radius_ball,...
                                                            [],[],...
                                                                type_value_func_computation,param);
                    else

                        value_func = KernelTCLValueFunc(number_of_points,time_horizon,...
                                                          type_vector_field,radius_ball,...
                                                            @GaussianKernel,kernel_parameter,...
                                                                type_value_func_computation,param);
                    end

                case 'KLdivAmbiguity'
                    fprintf('\nComputing value function (KL ambiguity)...\n');

                otherwise
                    not_implemented();
            end

        case 'LTI'
            switch struct_ambiguity_types.name
                case 'NoAmbiguity'

                    fprintf('\nComputing value function (without ambiguity)...\n');

                    number_of_points = param.number_of_points;
                    time_horizon = param.time_horizon;

                    value_func = NoAmbiguityLTIValueFunc(number_of_points,...
                        time_horizon,type_vector_field,param);
                case 'KernelAmbiguity'

                    fprintf('\nComputing value function (kernel ambiguity)...\n');

                    number_of_points = param.number_of_points;
                    time_horizon = param.time_horizon;

                    radius_ball = struct_ambiguity_types.radius_ball;
                    %kernel_parameter =
                    %struct_ambiguity_types.kernel_parameter;  THIS MAY NOT
                    %                                               BE NECESSARY
                    type_value_func_computation = struct_ambiguity_types.type_value_func_computation;

                    if strcmp(type_value_func_computation,'KME')
                        value_func = KernelLTIValueFunc(number_of_points,time_horizon,...
                                                          type_vector_field,radius_ball,...
                                                            [],[],...
                                                                type_value_func_computation,param);
                    end

                otherwise
                    not_implemented();
            end
% 
%                 case 'MomentAmbiguity'
%                     fprintf('\nComputing value function (moment ambiguity)...\n')
% 
%                     number_of_points = param.number_of_points;
%                     time_horizon = param.time_horizon;
%                     radius_mean = struct_ambiguity_types.radius_mean;
%                     radius_variance = struct_ambiguity_types.radius_variance;
% 
% 
%                     value_func = MomentTCLValueFunc(number_of_points,...
%                                         time_horizon,type_vector_field,...
%                                             radius_mean,radius_variance,...
%                                                     param);
%                 case 'KernelAmbiguity'
%                     fprintf('\nComputing value function (kernel ambiguity)...\n');
% 
%                     number_of_points = param.number_of_points;
%                     time_horizon = param.time_horizon;
% 
%                     radius_ball = struct_ambiguity_types.radius_ball;
%                     kernel_parameter = struct_ambiguity_types.kernel_parameter;
%                     type_value_func_computation = struct_ambiguity_types.type_value_func_computation;
% 
%                     if strcmp(type_value_func_computation,'KME')
%                         value_func = KernelTCLValueFunc(number_of_points,time_horizon,...
%                                                           type_vector_field,radius_ball,...
%                                                             [],[],...
%                                                                 type_value_func_computation,param);
%                     else
% 
%                         value_func = KernelTCLValueFunc(number_of_points,time_horizon,...
%                                                           type_vector_field,radius_ball,...
%                                                             @GaussianKernel,kernel_parameter,...
%                                                                 type_value_func_computation,param);
%                     end
% 
%                 case 'KLdivAmbiguity'
%                     fprintf('\nComputing value function (KL ambiguity)...\n');

%                 otherwise
%                     not_implemented();
%             end


        otherwise
            not_implemented();
    end

    value_func = value_func.get_index_safety(state_partition.get_values);
    value_func = value_func.backward_iteration(state_partition,...
                                                input_partition,...
                                                    param.outer_loop_info);

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