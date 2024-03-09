function out = run_LQR_example(LTI_system,obj_function,param)

add_function_paths(pwd)

% Defining system's matrices
A = LTI_system.A;
B = LTI_system.B;
initial_condition = LTI_system.initial_condition;

% Time horizon
time_horizon = param.time_horizon;

% Objective function
Q = obj_function.Q; R = obj_function.R;


% Parameters defining the ambiguity, etc
number_of_MC_simulations = param.number_of_MC_simulations;
mean_add_state_noise = param.mean_noise;
number_of_points = param.number_of_points;
number_of_input_partition = param.number_of_input_partition; 

radius_ball = param.radius_ball;
kernel_parameter = param.kernel_parameter;
type_value_func_computation = param.type_value_func_computation;
number_of_points_KME = param.number_of_points_KME;
regulariser_param = param.regulariser_param;
eta_param = param.eta_param;


chol_cov = chol(param.covaraiance_state_noise)';

discrete_sys = ss(A,B,eye(2),zeros(2,1),-1);

[gain_lqr,~,~] = lqr(discrete_sys,Q,R,[]);

LTI_vector_field = VectorFieldLTI(initial_condition,[],A,B);
out.LTI_sys = LTI_vector_field;

out.avg_obj_func_LQR = MC_compute_objective_func(LTI_vector_field,gain_lqr,...
                                   Q,R,mean_add_state_noise,chol_cov,...
                                        time_horizon,number_of_MC_simulations);


%%%%%%%%%%%%%%%% Generating trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

out.LQR_closed_loop_trajectory = generate_trajectory(LTI_vector_field,gain_lqr,mean_add_state_noise,chol_cov,time_horizon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        

internal_param = compute_transition('LTI',...
                                number_of_points,number_of_input_partition,number_of_MC_simulations); % Computing the transition probability

internal_param.time_horizon = int16(time_horizon);
internal_param.outer_loop_info = [];
internal_param.Q = Q;

state_partition = TwoDimStatePartition(internal_param.grid,number_of_points,internal_param.safe_set,'LTI');
input_partition = linspace(-5,5,number_of_input_partition)';

out.state_partition = state_partition;
out.input_partition = input_partition;

struct_ambiguity_types.name = 'NoAmbiguity';

time_kernel = tic;
value_func_no_ambiguity = main_value_function_iteration(state_partition,...
    input_partition,'LTI',struct_ambiguity_types,...
    exist('value_func_moment','var'),internal_param);

value_func_no_ambiguity.time = toc(time_kernel);

out.value_func_no_ambiguity = value_func_no_ambiguity;

%%%%%%%%%%%%%%%% Generating trajectories no-ambiguity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

out.No_ambiguity_closed_loop_trajectory = generate_trajectory_DP(state_partition,input_partition,LTI_vector_field,...
                                    value_func_no_ambiguity.opt_input(:,1),mean_add_state_noise,chol_cov,time_horizon);


%%%%%%%%%%%%%%% Kernel parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struct_ambiguity_types.name = 'KernelAmbiguity';
struct_ambiguity_types.radius_ball = radius_ball;
struct_ambiguity_types.kernel_parameter = kernel_parameter;
struct_ambiguity_types.type_value_func_computation = type_value_func_computation;
struct_ambiguity_types.number_of_samples_KME = number_of_points_KME;
struct_ambiguity_types.regulariser_param = regulariser_param;
struct_ambiguity_types.eta_param = eta_param;

internal_param.kernel_parameter = kernel_parameter;
internal_param.regulariser_param = regulariser_param;


[chol_fac,kernel_func,data_KME] = LTI_get_data_KME(state_partition,input_partition,...
    number_of_points_KME,internal_param);


internal_param.chol_fac = chol_fac;
internal_param.number_of_points_KME = number_of_points_KME;
internal_param.kernel_func = kernel_func;
internal_param.data_KME = data_KME;
internal_param.eta_param = eta_param;

outer_loop_info.string_ambiguity{1}.name = 'KernelAmbiguity';
outer_loop_info.type_vector_field = 'LTI';
outer_loop_info.total_iteration = 1;
outer_loop_info.current_iteration = 1;
outer_loop_info.time_iteration = 0;

internal_param.outer_loop_info = outer_loop_info;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time_kernel = tic;
value_func_kernel = main_value_function_iteration(state_partition,...
    input_partition,'LTI',struct_ambiguity_types,...
    exist('value_func_moment','var'),internal_param);

value_func_kernel.time = toc(time_kernel);

out.value_func_kernel = value_func_kernel;

%%%%%%%%%%%%%%%%%%%% Generating trajectory with Kernel policy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


out.kernel_closed_loop_trajectory = generate_trajectory_DP(state_partition,input_partition,LTI_vector_field,...
                                    value_func_kernel.opt_input(:,1),mean_add_state_noise,chol_cov,time_horizon);

end


%%%%%%%%%%%%%%%%%%%%% ADDITIONAL FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = generate_trajectory(LTI_vector_field_obj,controller,mean,chol_cov,time_horizon)
    
    state_trajectory = zeros(time_horizon+1,2);
    noise = zeros(time_horizon,2);
    control_action = zeros(time_horizon,1);

    state_trajectory(1,:) = LTI_vector_field_obj.get_values.current_state;

    for i = 1:time_horizon
        noise(i,:) = generate_noise_LTI(mean,chol_cov);
        control_action(i,:) = -controller*LTI_vector_field_obj.current_state;

        LTI_vector_field_obj.input = control_action(i,:);
        LTI_vector_field_obj.noise = noise(i,:)';

        state_trajectory(i+1,:) = LTI_vector_field_obj.iterate_dynamics()';

        LTI_vector_field_obj.current_state = state_trajectory(i+1,:)';
    end

    out.state = state_trajectory;
    out.noise = noise;
    out.control_action = control_action;

end

function out = generate_trajectory_DP(state_partition,input_partition,...
                                            LTI_vector_field_obj,controller,mean,chol_cov,time_horizon)
    
    state_trajectory = zeros(time_horizon+1,2);
    noise = zeros(time_horizon,2);
    control_action = zeros(time_horizon,1);

    state_trajectory(1,:) = LTI_vector_field_obj.get_values.current_state;

    for i = 1:time_horizon
        current_state_index = state_partition.get_element_partition(LTI_vector_field_obj.current_state).index_grid;
        noise(i,:) = generate_noise_LTI(mean,chol_cov);
        control_action(i,:) = input_partition(controller(current_state_index));

        LTI_vector_field_obj.input = control_action(i,:);
        LTI_vector_field_obj.noise = noise(i,:)';

        state_trajectory(i+1,:) = LTI_vector_field_obj.iterate_dynamics()';

        LTI_vector_field_obj.current_state = state_trajectory(i+1,:)';
    end

    out.state = state_trajectory;
    out.noise = noise;
    out.control_action = control_action;

end

function obj_func = compute_objective_func(Q,R,states_trajectories,control_actions)
    
    time_horizon = size(states_trajectories,1)-1;
    
    obj_func = 0;

    for i =1:time_horizon
        obj_func = obj_func + (states_trajectories(i,:)*Q*states_trajectories(i,:)' ...
                                    + control_actions(i,:)*R*control_actions(i,:));
    end


end

function estimate_obj_func = MC_compute_objective_func(LTI_vector_field_obj,controller,...
                                                            Q,R,mean,chol_cov,time_horizon,number_of_MC_simulations)

sum_MC = 0;

parfor i =1:number_of_MC_simulations
    temp = generate_trajectory(LTI_vector_field_obj,controller,mean,chol_cov,time_horizon);

    sum_MC = sum_MC + compute_objective_func(Q,R,temp.state,temp.control_action)/number_of_MC_simulations;
end

estimate_obj_func = sum_MC;

end


