clear all
close all

addpath AuxiliaryFunctions/
path_project = project_path;

add_function_paths(path_project);

% Empty to enable printing information of the script
struct_no_ambiguity.name = [];
struct_kernel_ambiguity.name = [];
struct_moment_ambiguity.name = [];
struct_KL_div_ambiguity.name = [];

time_horizon = int16(9);

% Number of points between 18 and 24 degree 
number_of_points = 150;
number_of_MC_simulation = 1000; % Set this to be a vector if want to run 
                              % the test with different value for this
                              % parameter
                              
number_of_samples = 200; % This is the number of samples for the kernel 
                         % ambiguity set. It may different from the number
                         % of samples used to estimate the transition
                         % probability

radius_ball = 0.1; % Define the radius of distance-based ambiguity sets.
                       % If this is a vector, it will run several
                       % simulations, one for each entry of the vector.

% The next two lines define values for the moment ambiguity set. If any of
% these are a vector, different simulations (once for each entry) are run.
radius_mean = [];
radius_variance = [];

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

% In addition to not passing as parameter in the function below, you should
% also comment here to ommit any ambiguity set
struct_no_ambiguity.name = 'NoAmbiguity';
struct_kernel_ambiguity.name = 'KernelAmbiguity';
struct_kernel_ambiguity.kernel_parameter = 20;
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

    TCL_results_path = [TCL_results_path;results_TCL];
    
    % If there exists distance-based ambiguity sets
    if ~isempty(radius_ball)
        struct_kernel_ambiguity.number_of_samples = number_of_samples;

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


paths_string = sprintf('./Results/results_%s/TCL/paths_%s.mat',...
                        char(java.net.InetAddress.getLocalHost.getHostName)...
                                ,TCL_results_path(end).file_name);

save(paths_string);

remove_function_paths(path_project);







