function out = TCL(time_horizon,struct_ambiguity_types,...
                        outer_loop_info,input_param)

param = input_param;

param.outer_loop_info = outer_loop_info;
param.outer_loop_info.string_ambiguity = struct_ambiguity_types; % THIS IS BAD PRACTICE SINCE I AM ALSO SHARING AMBIGUITY PARAMETERS
param.time_horizon = int8(time_horizon);

state_partition = OneDimStatePartition(param.grid,param.number_of_points,...
                                           param.safe_set,'TCL'); % Generate the partition 
                                                                  % of the state space
input_partition = generate_input_partition([],'TCL'); % Generate a vector with 
                                                      % all possible combinations of inputs

param_save.path_project = outer_loop_info.path_project;

number_of_ambiguity_sets = length(struct_ambiguity_types);

for i = 1:number_of_ambiguity_sets

    switch struct_ambiguity_types{i}.name

        case 'NoAmbiguity'
            TotalTime = tic;
            ValueFuncNoAmbiguity = MainValueFunctionIteration(state_partition,input_partition,'TCL',struct_ambiguity_types{i},exist('ValueFuncNoAmbiguity','var'),param);
            ValueFuncNoAmbiguity.time = toc(TotalTime);
            
            
        case 'MomentAmbiguity'
            time_moment_based = tic;
            value_func_moment = main_value_function_iteration(state_partition,input_partition,'TCL',struct_ambiguity_types{i},exist('value_func_moment','var'),param);
            value_func_moment.time = toc(time_moment_based);
            
            param_save.radius_mean = struct_ambiguity_types{i}.radius_mean;
            param_save.radius_variance = struct_ambiguity_types{i}.radius_variance;
             
        case 'WassersteinAmbiguity'
            TotalTime2 = tic;
            ValueFuncWasserstein = MainValueFunctionIteration(state_partition,input_partition,'TCL',struct_ambiguity_types{i},exist('ValueFuncWasserstein','var'),param);
            ValueFuncWasserstein.time = toc(TotalTime2);
            
            param_save.ep = struct_ambiguity_types{i}.ep;
                   
            
        case 'KLdivAmbiguity'
            TotalTime3 = tic;
            ValueFuncKL = MainValueFunctionIteration(state_partition,input_partition,'TCL',struct_ambiguity_types{i},exist('ValueFuncKL','var'),param);
            ValueFuncKL.time = toc(TotalTime3);
            
            param_save.ep = struct_ambiguity_types{i}.ep;
            
        case 'KernelAmbiguity'  
            TotalTime4 = tic;
            ValueFuncKernel = MainValueFunctionIteration(state_partition,input_partition,'TCL',struct_ambiguity_types{i},exist('ValueFuncKernel','var'),param);
            ValueFuncKernel.time = toc(TotalTime4);
            
            param_save.ep = struct_ambiguity_types{i}.ep;
           
        otherwise
            warning('%s has not been implemented. Jumping to the next string',struct_ambiguity_types{i});
    end   
end

file_name = get_date_save_file(time_horizon,param.number_of_points,...
                            param.number_of_MC_simulations,'TCL',...
                                param_save); % Getting the name of file based on the current date and time

save(file_name.full_path); % saving the results in the path specified by FILE

out = file_name;

end

