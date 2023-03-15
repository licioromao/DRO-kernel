clear all
addpath('./VectorFields');
addpath('./PartitionClassesAndFunctions/')
addpath('./EstimatingTransProb/')
addpath('./ValueFunctionComputation/')

% Standard LQR control
close all
A = [1,0;1,1]; B = [1;1]; initial_condition = 10*[1;2]; time_horizon = 1000;
Q = eye(2); R = 1;
number_of_MC_simulations = 1000;

covaraiance_state_noise = randn(2,1);
mean_add_state_noise = [0;0]; covaraiance_state_noise = covaraiance_state_noise*covaraiance_state_noise' + 0.1*eye(2);
chol_cov = chol(covaraiance_state_noise)';

%%%%%%%%%%%%%%%% Testing noise function %%%%%%%%%%%%%%%%

% m = 100000;
% temp_noise = zeros(2,m);
% for i=1:m
%     temp_noise(:,i) = generate_noise_LTI(mean_add_state_noise,chol_cov);
% end
% 
% estimate_mean = sum(temp_noise,2)/m;
% estimate_cov = (temp_noise-estimate_mean)*(temp_noise-estimate_mean)'/m;
% 
% [estimate_cov,covaraiance_state_noise]
% 
% [estimate_mean,mean_add_state_noise]

%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%

discrete_sys = ss(A,B,eye(2),zeros(2,1),-1);

[gain_lqr,solution_ricatti,closed_loop_poles] = lqr(discrete_sys,Q,R,[]);

LTI_vector_field = VectorFieldLTI(initial_condition,[],A,B);

MC_compute_objective_func(LTI_vector_field,gain_lqr,...
                                   Q,R,mean_add_state_noise,chol_cov,...
                                        time_horizon,number_of_MC_simulations)


%%%%%%%%%%%%%%%% Generating trajectories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

temp = generate_trajectory(LTI_vector_field,gain_lqr,mean_add_state_noise,chol_cov,time_horizon);
temp2 = generate_trajectory(LTI_vector_field,gain_lqr,mean_add_state_noise,chol_cov,time_horizon);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

number_of_points = [3,3];
number_of_MC_simulation = 10;

param = compute_transition('LTI',...
                                number_of_points,number_of_MC_simulation); % Computing the transition probability

param.time_horizon = time_horizon;

state_partition = TwoDimStatePartition(param.grid,number_of_points,param.safe_set,'LTI');
input_partition = linspace(-5,5,20)';

time_no_ambiguity = tic;
value_func_no_ambiguity = main_value_function_iteration(state_partition,...
    input_partition,'LTI',[],...
    exist('value_func_moment','var'),param);

value_func_no_ambiguity.time = toc(time_no_ambiguity);


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

function obj_func = compute_objective_func(Q,R,states_trajectories,control_actions)
    
    time_horizon = size(states_trajectories,1)-1;
    
    obj_func = 0;

    for i =1:time_horizon
        obj_func = obj_func + (states_trajectories(i,:)*Q*states_trajectories(i,:)' ...
                                    + control_actions(i,:)*R*control_actions(i,:))/time_horizon;
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


