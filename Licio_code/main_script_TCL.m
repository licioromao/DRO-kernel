clear all

time_horizon = 9;
number_of_grid_points = 100;
number_of_MC_transitions = 1000;
number_of_samples_KME = 150;

radius_ball = [0,0.01,0.05];

kernel_parameters = [1,20,1];
%                       1,1,1;
%                         10,1,1;...
%                          100,1,1;...
%                            10,0.2,1;...
%                              10,2,1;...
%                                10,20,1;...
%                                  10,200,1];

number_of_simulations = size(kernel_parameters,1);
results = cell(number_of_simulations,1);



for i=1:number_of_simulations
    results{i} = testTCLFunc(time_horizon,number_of_grid_points,number_of_MC_transitions,...
                                number_of_samples_KME,radius_ball,[],[],pwd,kernel_parameters(i,:)); 
end

save test_TCL_2.mat