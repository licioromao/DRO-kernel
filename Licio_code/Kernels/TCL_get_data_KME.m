function [chol_fac,kernel_func,data_KME] = ...
                    TCL_get_data_KME(obj_state_partition,input_partition...
                                            ,number_of_points_KME,param)

grid_no_inputs = obj_state_partition.get_values.partition.X;
number_input_partition = size(input_partition,1);

grid_with_inputs = sort(repmat(grid_no_inputs,number_input_partition,1));
grid_with_inputs = [grid_with_inputs,repmat(input_partition,...
                        size(grid_no_inputs,1),1)]; % Forming state-input partition

number_state_input_grid = size(grid_with_inputs,1);
input_data = sort(randi(number_state_input_grid,[number_of_points_KME,1]));

data = grid_with_inputs(input_data,:);

y = zeros(number_of_points_KME,1);

for i = 1:number_of_points_KME
    temp = VectorFieldTCL(data(i,1),data(i,2),param);
    noise = generate_noise(param,'TCL');
    temp.noise = noise;
    
    y(i) = temp.iterate_dynamics();
end

data = [data,y]; % Collecting state-input-next_state data

gram_matrix = GaussianKernel(data(:,1),data(:,1),param.kernel_parameter,...
                                    [],'KME',data(:,3));

chol_fac = chol(gram_matrix + param.regulariser_param*...
                    number_of_points_KME*eye(number_of_points_KME)); % Cholesk factorisation of Gram matrix
kernel_func = @(state,input)compute_kernel_regression(data,state,input,param.kernel_parameter); 
                                % Function to compute embedding for a given
                                % state-input pair

data_KME = data;
end

function out = compute_kernel_regression(data,state,input,kernel_parameter)

number_of_points = size(data,1);
out = zeros(number_of_points,1);

for i=1:number_of_points 
    out(i) = GaussianKernel(data(i,1),state,kernel_parameter,[],'KME_2',[data(i,3),input]);
end

end