clear all

time_horizon = 18;
number_of_grid_points = 35;
number_of_MC_transitions = 1000;
number_of_samples_KME = 7000;

radius_ball = [0,0.01,0.05];
% 
kernel_parameters = [100,200,1;
                        100,200,0.97];
%                             100,200,0.9;
%                                 100,200,0.6;
%                                     100,200,0.3];
%                         100,2000,1;
%                             100,20000,1;
%                                 100,200000,1];
%                         100,500,1;
%                             100,1000,1;
%                             1000,1000,1;
%                             10000,1000,1];
%                                 300,200,1;
%                                    300,500,1;
%                                      300,1000,1];

radius_mean = [0.2,0.3];
radius_variance = [1.5,3.5];

% [10,200,1;
%                         10,500,1;
%                             10,1000,1;
%                                 0.10,200,1;
%                                   .10,500,1;
%                                      .10,1000,1];
%                        1,1,1;
%                          10,1,1;...
%                           100,1,1;...
%                             10,0.2,1;...
%                               10,2,1;...
%                                 10,20,1;...
%                                   10,200,1];

number_of_simulations = size(kernel_parameters,1);
results = cell(number_of_simulations,1);



for i=1:number_of_simulations
%     results{i} = testTCLFunc(time_horizon,number_of_grid_points,number_of_MC_transitions,...
%                                 number_of_samples_KME,radius_ball,[],[],pwd,kernel_parameters(i,:)); 
   results{i} = testTCLFunc(time_horizon,number_of_grid_points,number_of_MC_transitions,...
                               number_of_samples_KME,radius_ball,radius_mean,radius_variance,pwd,kernel_parameters(i,:)); 
end

save CDC_2.mat