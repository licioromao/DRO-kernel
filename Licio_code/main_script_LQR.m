
%% Parameters of the model

number_of_MC_simulations = 300;
mean_add_state_noise = [0;0];
covaraiance_state_noise = randn(2,1);
covaraiance_state_noise = covaraiance_state_noise*covaraiance_state_noise' + eye(2);

number_of_points = [25,25];
number_of_input_partition = 50;

radius_ball = 0.01;
kernel_parameter = 0.001;
type_value_func_computation = 'KME';
number_of_points_KME = 800;
regulariser_param = 0.2;
eta_param = 1;

time_horizon = 15;


% Number of discretisation points
param.number_of_MC_simulations = number_of_MC_simulations;
param.mean_noise = mean_add_state_noise;
param.covaraiance_state_noise = covaraiance_state_noise;
param.number_of_points = number_of_points;
param.number_of_input_partition = number_of_input_partition; 

% Kernel parameters
param.radius_ball = radius_ball;
param.kernel_parameter = kernel_parameter;
param.type_value_func_computation = type_value_func_computation;
param.number_of_points_KME = number_of_points_KME;
param.regulariser_param =  regulariser_param;
param.eta_param = eta_param;

% Time-horizon
param.time_horizon = time_horizon;

%% LTI system

initial_condition = [5;3];

LTI_sys.A = [1,0;1,1]; LTI_sys.B = [1;1]; LTI_sys.initial_condition = initial_condition;

%% Objective function

obj_func.Q = 10*eye(2); obj_func.R = 1;


aa = run_LQR_example(LTI_sys,obj_func,param);

[aa.LQR_closed_loop_trajectory.state(1:end-1,:),aa.LQR_closed_loop_trajectory.control_action]

[aa.kernel_closed_loop_trajectory.state(1:end-1,:),aa.kernel_closed_loop_trajectory.control_action]

[aa.No_ambiguity_closed_loop_trajectory.state(1:end-1,:),aa.No_ambiguity_closed_loop_trajectory.control_action]


