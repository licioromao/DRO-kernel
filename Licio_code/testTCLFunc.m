function out = testTCLFunc(horizon,number_of_points,number_of_MC_simulation...
                            ,number_of_samples_KME,radius_ball,radius_mean,...
                                 radius_variance,path_project,kernel_parameters)

time_horizon = int16(horizon); 

% Empty to enable printing information of the script
struct_no_ambiguity.name = [];
struct_kernel_ambiguity.name = [];
struct_moment_ambiguity.name = [];
%struct_KL_div_ambiguity.name = [];


if ~isempty(number_of_MC_simulation)
    number_of_sumulations_TCL = size(number_of_MC_simulation,2);
    total_iterations = number_of_sumulations_TCL;
else
    number_of_sumulations_TCL = 0;
    total_iterations = 1;
end


if ~isempty(radius_ball)
    number_of_radius_distance = size(radius_ball,2);
    total_iterations = total_iterations*number_of_radius_distance;
else
    number_of_radius_distance = 0;
end


if ~isempty(radius_mean)
    number_of_simulations_mean_moment = size(radius_mean,2);
    total_iterations = total_iterations*number_of_simulations_mean_moment;
else
    number_of_simulations_mean_moment = 0;
end


if ~isempty(radius_variance)
    number_of_simulations_variance_moment = size(radius_variance,2);
    total_iterations = total_iterations*number_of_simulations_variance_moment;
else
    number_of_simulations_variance_moment = 0;
end

TCL_results_path = cell(total_iterations,1);

% In addition to not passing as parameter in the function below, you should
% also comment here to ommit any ambiguity set
struct_no_ambiguity.name = 'NoAmbiguity';
struct_kernel_ambiguity.name = 'KernelAmbiguity';
struct_kernel_ambiguity.kernel_parameter = kernel_parameters(1);
struct_kernel_ambiguity.regulariser_param = kernel_parameters(2);
struct_kernel_ambiguity.eta_param = kernel_parameters(3);

struct_moment_ambiguity.name = 'MomentAmbiguity';
%struct_KL_div_ambiguity.name = 'KLdivAmbiguity';


time_per_iteration = 2; % A guess on the first iteration


outer_loop_info.type_vector_field = 'TCL';
outer_loop_info.path_project = path_project;
outer_loop_info.total_iteration = total_iterations;
outer_loop_info.time_iteration = time_per_iteration;

TCL_results_path = [];

for i1=1:number_of_sumulations_TCL
    param = compute_transition(outer_loop_info.type_vector_field,...
                                number_of_points,number_of_MC_simulation(i1));

    index = remaining_iterations(1,[i1,number_of_sumulations_TCL],...
                                    total_iterations/number_of_sumulations_TCL,[]);

    outer_loop_info.current_iteration = index;

    results_TCL = TCL(time_horizon,{struct_no_ambiguity},...
                        outer_loop_info,param); % Solve the DP iteration without ambiguity set
    
    pause(2)
    TCL_results_path = [TCL_results_path;results_TCL];
    
    % If there exists distance-based ambiguity sets
    if ~isempty(radius_ball)
        struct_kernel_ambiguity.number_of_samples_KME = number_of_samples_KME;


        for i2=1:number_of_radius_distance
            % Parameters of the kernel ambiguity set
            initial_time_per_iteration = tic;
            struct_kernel_ambiguity.radius_ball = radius_ball(i2);
            index = remaining_iterations(2,[[i1;i2],...
                                [number_of_sumulations_TCL;number_of_radius_distance]]...
                                ,number_of_simulations_mean_moment...
                                    *number_of_simulations_variance_moment,[]);
            outer_loop_info.current_iteration = index;

            results_TCL = TCL(time_horizon,{struct_kernel_ambiguity},...
                                outer_loop_info,param); % Solving DP with Kernel Ambiguity set
            TCL_results_path = [TCL_results_path;results_TCL];
            pause(2)

            time_per_iteration = toc(initial_time_per_iteration); % Estimating how long time it took the last iteration
        end
    end

    if ~isempty(radius_mean) && ~isempty(radius_variance)

        for i3=1:number_of_simulations_mean_moment
            for i4=1:number_of_simulations_variance_moment
                

                % Parameters of the moment ambiguity set
                struct_moment_ambiguity.radius_mean = radius_mean(i3);
                struct_moment_ambiguity.radius_variance = radius_variance(i4);

                % Parameters of the kernel ambiguity set
                initial_time_per_iteration = tic;

                index = remaining_iterations(2,[[i3;i4],...
                        [number_of_simulations_mean_moment;...
                            number_of_simulations_variance_moment]],1,[]);

                outer_loop_info.CurrentIteration = index;

                results_TCL = TCL(time_horizon,{struct_moment_ambiguity},...
                                    outer_loop_info,param); % Solving DP problem with Moment ambiguity set
                
                TCL_results_path = [TCL_results_path;results_TCL];

                time_per_iteration = toc(initial_time_per_iteration); % Estimating how long time it took the last iteration

            end
        end
    end
end

paths_string = sprintf('/Results/results_%s/TCL/paths_%s.mat',...
                        char(java.net.InetAddress.getLocalHost.getHostName)...
                                ,TCL_results_path(end).file_name);

out.paths_internal = strcat(path_project,paths_string);
out.paths_results = TCL_results_path;

save(out.paths_internal)

end