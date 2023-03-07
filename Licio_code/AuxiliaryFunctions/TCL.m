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
            time_no_ambiguity = tic;
            value_func_no_ambiguity = main_value_function_iteration(state_partition,...
                                            input_partition,'TCL',struct_ambiguity_types{i},...
                                                exist('value_func_moment','var'),param);

            value_func_no_ambiguity.time = toc(time_no_ambiguity);
            
        case 'MomentAmbiguity'
            time_moment_based = tic;
            value_func_moment = main_value_function_iteration(state_partition,...
                                        input_partition,'TCL',struct_ambiguity_types{i},...
                                                exist('value_func_moment','var'),param);

            value_func_moment.time = toc(time_moment_based);
            
            param_save.radius_mean = struct_ambiguity_types{i}.radius_mean;
            param_save.radius_variance = struct_ambiguity_types{i}.radius_variance;
            
        case 'KernelAmbiguity'
            
            % The parameters below are necessary to perform kernel mean
            % embedding
            param.regulariser_param = [];
            param.eta_param = [];
            param.kernel_parameter = [];

            param.chol_fac = [];
            param.number_of_points_KME = [];
            param.kernel_func = [];

%             % Conservative computation of the kernel
%             time_kernel_based = tic;
%             struct_ambiguity_types{i}.type_value_func_computation = 'Conservative';
% 
%             value_func_kernel_conservative = main_value_function_iteration(state_partition,...
%                 input_partition,'TCL',struct_ambiguity_types{i},...
%                 exist('value_func_kernel_conservative','var'),param);
% 
%             value_func_kernel_conservative.time = toc(time_kernel_based);

%             % Our second idea to compute the value function under kernel
%             % ambiguity
% 
%             time_kernel_based = tic;
%             struct_ambiguity_types{i}.type_value_func_computation = 'QP';
% 
%             value_func_kernel_QP = main_value_function_iteration(state_partition,...
%                 input_partition,'TCL',struct_ambiguity_types{i},...
%                 exist('value_func_kernel_QP','var'),param);
% 
%             value_func_kernel_QP.time = toc(time_kernel_based);
%              
%               
%               % The computation below solves a regression problem for the
%               % kernel
% 
%             time_kernel_based = tic;
%             struct_ambiguity_types{i}.type_value_func_computation = 'Matrix';
% 
%             value_func_kernel_matrix = main_value_function_iteration(state_partition,...
%                 input_partition,'TCL',struct_ambiguity_types{i},...
%                 exist('value_func_kernel_matrix','var'),param);
% 
%             value_func_kernel_matrix.time = toc(time_kernel_based);

            % Our last attempt using the kernel mean embedding

            time_kernel_based = tic;
            
            number_of_samples_KME = struct_ambiguity_types{i}.number_of_samples_KME;
            
            struct_ambiguity_types{i}.type_value_func_computation = 'KME';
            param.regulariser_param = struct_ambiguity_types{i}.regulariser_param;
            param.eta_param = struct_ambiguity_types{i}.eta_param;
            param.kernel_parameter = struct_ambiguity_types{i}.kernel_parameter;

            [chol_fac,kernel_func,data_KME] = TCL_get_data_KME(state_partition,input_partition,...
                                    number_of_samples_KME,param);
            
           
            param.chol_fac = chol_fac;
            param.number_of_samples_KME = number_of_samples_KME;
            param.kernel_func = kernel_func;
            param.data_KME = data_KME;

            value_func_kernel_KME = main_value_function_iteration(state_partition,...
                                    input_partition,'TCL',struct_ambiguity_types{i},...
                                            exist('value_func_kernel_KME','var'),param);

            value_func_kernel_KME.time = toc(time_kernel_based);


            param_save.radius_ball = struct_ambiguity_types{i}.radius_ball;
            param_save.kernel_parameter = struct_ambiguity_types{i}.kernel_parameter;

        case 'KLdivAmbiguity'
            not_implemented();

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

