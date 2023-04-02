function [] = plot_TCL_results(output_TCL_func,project_path,number_of_MC_simulations,save,filename)

path_to = '/Users/licioromao/OneDrive - Nexus365/Postdoc/Papers/Ashish_collaboration/CDC-paper/ACC23_AshishLicio_KernelSafety/Figures/';

add_function_paths(project_path)

path_results = output_TCL_func.paths_results;
number_of_TCL_experiments = size(path_results,1);

% Retriving information from saved files
variables = {'state_partition','input_partition','param','value_func*'};

for i =1:number_of_TCL_experiments
    if exist('value_func_kernel_KME','var')
        new_var = change_var_name(value_func_kernel_KME,'KernelAmbiguity');
        eval(new_var);
    end

    if exist('value_func_moment','var')
        new_var = change_var_name(value_func_moment,'MomentAmbiguity');
        eval(new_var);
    end

    load(path_results(i).full_path,variables{:});
    whos
end

if exist('value_func_kernel_KME','var')
    new_var = change_var_name(value_func_kernel_KME,'KernelAmbiguity');
    eval(new_var);
    clear value_func_kernel_KME
end

if exist('value_func_moment','var')
    new_var = change_var_name(value_func_moment,'MomentAmbiguity');
    eval(new_var);
    clear value_func_moment
end

% Plotting results from the kernel ambiguity
results_kernel = whos('value_func_kernel_KME_*');
number_of_KME = size(results_kernel,1);

if number_of_KME > 0
    figure
    hold on
    name_legend_KME = [];

    for i =1:number_of_KME
        plot(state_partition.partition.grid_x,...
            eval(results_kernel(i).name).value_function(:,1),'LineWidth',2);

        name_legend_KME = [name_legend_KME,{sprintf('\\epsilon = %.4f', ...
            eval(results_kernel(i).name).radius_ball)}];
    end

    if ~isfield(param,'kernel_parameter')
        param.kernel_parameter = eval(results_kernel(i).name).param.kernel_parameter;
        param.regulariser_param = eval(results_kernel(i).name).param.regulariser_param;
        param.eta_param = eval(results_kernel(i).name).param.eta_param;
    end

    legend(name_legend_KME)
    box on
    grid on
    text_title = {'Safety probability. Kernel Ambiguity',sprintf('Kernel parameter: %.4f,  Regulariser: %.4f,  Eta:%.4f',...
                                                        param.kernel_parameter,param.regulariser_param,param.eta_param)};
    if save
        title(text_title,'Interpreter','latex')
        individual_plot = gcf;s
        temp = strcat(path_to,filename,'-all','.eps');
        exportgraphics(individual_plot,temp)
    end
end

% Plotting results from the moment ambiguity
results_moment = whos('value_func_moment*');
number_of_moment = size(results_moment);

if number_of_moment > 0
    figure
    hold on
    name_legend_moment = [];

    for i=1:number_of_moment
        plot(state_partition.partition.grid_x,...
            eval(results_moment(i).name).value_function(:,1),'LineWidth',2);

        name_legend_moment = [name_legend_moment,{sprintf('\\mu = %.4f, \\Sigma = %.4f', ...
            eval(results_moment(i).name).radius_mean, ...
            eval(results_moment(i).name).radius_variance)}];
    end
    legend(name_legend_moment)
    box on
    grid on
    text_title = {'Safety probability -- Moment Ambiguity sets'};

    title(text_title,'Interpreter','latex')

    if save
        moment_plot = gcf;
        temp = strcat(path_to,filename,'-moment','.eps');
        exportgraphics(moment_plot,temp);
    end
    
end

% Performing Monte Carlo simulation of the designed policy, if last
% parameter is not empty

if ~isempty(number_of_MC_simulations)
    
    grid_state = state_partition.partition.grid_x;
    number_of_points = length(grid_state);
    
    empirical_value_function_KME = zeros(number_of_points,number_of_KME);
    % Iterating over the kernel value function
    for i=1:number_of_KME
        opt_input = eval(results_kernel(i).name).opt_input;
        time_horizon = eval(results_kernel(i).name).time_horizon;
        param = eval(results_kernel(i).name).param;

        for j = 1:number_of_points
            current_state = grid_state(j);
            empirical_value_function_KME(j,i) = compute_empirical_probability(current_state,...
                                    state_partition,opt_input,time_horizon,...
                                        number_of_MC_simulations,param);
        end
        
        figure
        hold on
        plot(state_partition.partition.grid_x,...
            eval(results_kernel(i).name).value_function(:,1),'LineWidth',2);
        plot(state_partition.partition.grid_x,...
            empirical_value_function_KME(:,i),'k--','LineWidth',2);

        name_legend_KME = {sprintf('\\epsilon = %.4f',...
                                eval(results_kernel(i).name).radius_ball),'Empirical'};

        legend(name_legend_KME)
        box on
        grid on
        text_title = sprintf('Kernel parameter: %.4f,  Regulariser: %.4f,  Eta:%.4f',...
                                                        param.kernel_parameter,param.regulariser_param,param.eta_param);
        title({'Value function versus empirical value function',text_title},'Interpreter','latex');

        if save
            all_plot = gcf;
            temp = strcat(path_to,filename,'-individual-',sprintf('%d',i),'.eps');
            exportgraphics(all_plot,temp)
        end
    end
end

end

function out = change_var_name(value_func,ambiguity_type)

switch ambiguity_type
    case 'KernelAmbiguity'

        radius_ball = value_func.radius_ball;            
        temp = strcat('value_func_kernel_KME_radius_',sprintf('%.4f',radius_ball));
        temp(strfind(temp,'.')) = '_';
        out = [temp,'=','value_func_kernel_KME'];

    case 'MomentAmbiguity'

        radius_mean = value_func.radius_mean;
        radius_variance = value_func.radius_variance;

        temp = strcat('value_func_moment_mean_',sprintf('%.4f',radius_mean));
        temp = strcat(temp,'_variance_',sprintf('%.4f',radius_variance));
        temp(strfind(temp,'.')) = '_';
        out = [temp,'=','value_func_moment'];

    otherwise
        not_implemented();
end

end

function out = simulate_trajectory(current_state,state_partition,opt_input,...
                                    time_horizon,param)

safe_set = param.safe_set;
out = 1;

temp_parition = state_partition.get_element_partition(current_state);
u = opt_input(temp_parition.index);

TCL_vector_field = VectorFieldTCL(current_state,u,param);

if current_state < safe_set(1) || current_state > safe_set(2)
    out = 0;
    return;
else
    for k = 1:time_horizon
        TCL_vector_field.noise = random(param.w);
        current_state = TCL_vector_field.iterate_dynamics();

        if current_state < safe_set(1) || current_state > safe_set(2)
            out = 0;
            return;
        end

        temp_parition = state_partition.get_element_partition(current_state);
        u = opt_input(temp_parition.index);

        TCL_vector_field.set_values(current_state,u,param);
    end
end

end

function out = compute_empirical_probability(current_state,state_partition,...
                                                opt_input,time_horizon,...
                                                    number_of_MC_simulations,param)

out = 0;

for k = 1:number_of_MC_simulations
    out = out + simulate_trajectory(current_state,state_partition,opt_input...
                                            ,time_horizon,param);
end

out = out/number_of_MC_simulations;

end