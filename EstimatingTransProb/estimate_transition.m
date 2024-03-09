function out = estimate_transition(obj_state_partition,input_partition,param)

% This function iterates over all state-action pairs and computes the
% transition probability.
%
%
%   Input: obj_one_dim -- This is a OneDimStateParition object.
%          input_partition -- output of the generateInput partition function
%          param -- a structure with all parameters of the problem
%
%
%   Output: TBD


type_vector_field = obj_state_partition.type_vector_field;

grid_with_inputs = param.grid_with_inputs;

switch type_vector_field
    
    case 'TCL'
        
        out = estimate_prob_transition(obj_state_partition,input_partition,...
                                        grid_with_inputs,'1D',type_vector_field,param);

    case 'LTI'

        out = estimate_prob_transition(obj_state_partition,input_partition,...
                                        grid_with_inputs,'2D',type_vector_field,param);
        
%     case 'ChainInt'
%         
%         out = estimate_prob_transition(obj_state_partition,input_partition,...
%                                         grid_no_inputs,'2D',type_vector_field,param);
%     
%     case 'Fishery'
%         
%         out = estimate_prob_transition(obj_state_partition,input_partition,...
%                                         grid_no_inputs,'3D',type_vector_field,param);
%     
%     case 'CarPole'
%         
%         out = estimate_prob_transition(obj_state_partition,input_partition,...
%                                         grid_no_inputs,'4D',type_vector_field,param);
%         
%     case 'CarPoleNL'
%         
%         out = estimate_prob_transition(obj_state_partition,input_partition,...  
%                                         grid_no_inputs,'4D',type_vector_field,param);
        
    otherwise
        not_implemented();
end

end


% Computing the transition probabilities

function out = estimate_prob_transition(obj_state_partition,input_partition,...
                                           grid_with_inputs,dim_space,type_vector_field,param)
    
    number_parallel_cores = 12; % Adjust this parameter based on your PC.
                                % Used only for parallel computations.
                                % THIS COULD BE MADE PC SPECIFIC

    temp_state_partition = obj_state_partition.get_values.partition;
    
    switch dim_space
        case '1D'
            
            number_of_points = size(temp_state_partition.X,1); % Number of points in the first dimension
            number_of_inputs = length(input_partition); % Number of possible inputs
            
            number_of_MC_simulations = param.number_of_MC_simulations; % Number of Monte Carlo simulation for each transition
            
            total_iterations = number_of_points*number_of_inputs; % remaining iterations
            temp_values = cell(total_iterations,1);
            print_estimate_transition_prob(total_iterations,0,0);
            

            % Variables to estimate the remaining time
            sum_time = 0;
            sum_it = 0;
            
            for i = 1:number_of_points
                x = temp_state_partition.X(i); % Get the current state
                
                if number_of_MC_simulations <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                    
                    % Iterate over all possible inputs
                    intial_time_input = tic;
                    temp_prob = zeros(number_of_points,number_of_inputs);
                    
                    for j =1:number_of_inputs
                        u = input_partition(j,:); % selecting a particular input
                        prob_xu = run_monte_carlo(x,u,obj_state_partition,...
                                                    type_vector_field,param); % empirical estimate of the transition probability
                        temp_prob(:,j) = prob_xu; % storing estimate
                    end
                    
                    index = remaining_iterations(1,[i,number_of_points],...
                                                    number_of_inputs,[]); % remaining iterates. It takes into account the iterates of the previous for loop
                    
                    for j=number_of_inputs-1:-1:0
                        temp_values{index - j} = temp_prob(:,number_of_inputs-j);
                    end
                    
                    final_time_input = toc(intial_time_input);
                    sum_time = sum_time + final_time_input;
                    sum_it = sum_it + 1;
                      
                    % Printing on the screen. THIS SHOULD BECOME A FUNCTION
                    if (number_of_points-1)/number_of_inputs < 100
                        if mod(i+j,5) == 0
                            print_estimate_transition_prob(total_iterations,...
                                                            index,sum_time/sum_it);
                            sum_time = 0;
                            sum_it = 0;
                        end
                    else
                        if mod(i+j,15) == 0
                            print_estimate_transition_prob(total_iterations,...
                                                            index,sum_time/sum_it);
                            sum_time = 0;
                            sum_it = 0;
                        end
                    end
                    
                else % if number of MC simulation is larger than 100, we will 
                     % leverage parallel computation
                    
                    % Iterate over all possible input combination
                    for j =1:number_of_inputs
                        u = input_partition(j); % selecting a particular 
                                                % allowable input
                        
                        initial_time_input = tic;
                        
                        prob_xu = run_monte_carlo_parallel(x,u,obj_state_partition,...
                                                            type_vector_field,number_parallel_cores,param); % empirical estimate of the transition probability using parallel computation
                        
                        index = remaining_iterations(2,[[i;j],...
                                        [number_of_points;number_of_inputs]],1,[]);
                        temp_values{index} = prob_xu;
                        
                        final_time_input = toc(initial_time_input);
                        sum_time = sum_time + final_time_input;
                        sum_it = sum_it + 1;
                        
                        % Printing on the screen. THIS SHOULD BE A FUNCTION
                        if (number_of_points-1)/number_of_inputs < 100
                            if mod(i+j,5) == 0
                                print_estimate_transition_prob(total_iterations,...
                                                              index,sum_time/sum_it);
                                sum_time = 0;
                                sum_it = 0;
                            end
                        else
                            if mod(i+j,15) == 0
                                print_estimate_transition_prob(total_iterations,...
                                                                index,sum_time/sum_it);
                                sum_time = 0;
                                sum_it = 0;
                            end
                        end
                    end
                end
            end
            
        % STOPPED HERE.    
        case '2D'
            
            Nx1 = size(temp_state_partition.X1,2); % Number of points in the 
                                                   % first dimension
            Nx2 = size(temp_state_partition.X1,1); % Number of points in the 
                                                   % second dimension
            
            number_of_inputs = length(input_partition); % Number of possible inputs
            
            number_of_MC_simulations = param.number_of_MC_simulations; % Number of Monte Carlo simulation for each transition
            
            number_of_points = Nx1*Nx2;
            temp_values = cell(Nx1*Nx2*number_of_inputs,1);
            
            total_iterations = Nx1*Nx2*number_of_inputs; % total number of 
                                                         % remaining iterations
            print_estimate_transition_prob(total_iterations,0,0);
            
            % Variable to estimate the amount of remaining time
            sum_time = 0;
            sum_it = 0;
            
            
            for i1 = 1:Nx1
                for i2 = 1:Nx2
                    
                    x = [temp_state_partition.X1(i2,i1);...
                              temp_state_partition.X2(i2,i1)]; % Get the current state
                    if number_of_MC_simulations <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                        
                        % Iterate over all possible input combination
                        tempTime = tic;
                        temp_prob = zeros(number_of_points,number_of_inputs);
                        for j =1:number_of_inputs % THIS SHOULD BE PARFOR
                            u = input_partition(j,:); % selecting a particular 
                                                      % allowable input
                            prob_xu = run_monte_carlo(x,u,obj_state_partition,...
                                   type_vector_field,param); % empirical estimate 
                                                             % of the transition probability
                            temp_prob(:,j) = prob_xu;
                        end
                        
                        index = remaining_iterations(2,[[i1;i2],[Nx1;Nx2]],...
                                                           number_of_inputs,[]);
                        
                        for j=number_of_inputs-1:-1:0
                            temp_values{index - j} = temp_prob(:,number_of_inputs-j);
                        end
                        
                        final_time_input = toc(tempTime);
                        sum_time = sum_time + final_time_input;
                        sum_it = sum_it + 1;
                        
                        % Printing results on the screen
                        if (number_of_points-1)/number_of_inputs < 100
                            if mod(i1+i2,5) == 0
                                print_estimate_transition_prob(total_iterations,index,sum_time/sum_it);
                                sum_time = 0;
                                sum_it = 0;
                            end
                        else
                            if mod(i1+i2,15) == 0
                                print_estimate_transition_prob(total_iterations,index,sum_time/sum_it);
                                sum_time = 0;
                                sum_it = 0;
                            end
                        end
                        
                    else % if number of MC simulation is larger than 100, we will leverage parallel computation
                        
                        % Iterate over all possible input combination
                        sum_time = 0;
                        sum_it = 0;
                        for j =1:number_of_inputs
                            u = input_partition(j,:); % selecting a particular allowable input
                            
                            tempTime1 = tic;
                            
                            prob_xu = run_monte_carlo_parallel(x,u,obj_state_partition,type_vector_field,12,param); % empirical estimate of the transition probability using parallel computation
                            
                            index = remaining_iterations(3,[[i1;i2;j],[Nx1;Nx2;number_of_inputs]],1,[]);
                            temp_values{index} = prob_xu;
                            
                            final_time_input = toc(tempTime1);
                            sum_time = sum_time + final_time_input;
                            sum_it = sum_it + 1;
                            
                            % Printing on the screen
                            if (number_of_points-1)/number_of_inputs < 100
                                if mod(i1+i2,5) == 0
                                    print_estimate_transition_prob(total_iterations,index,sum_time/sum_it);
                                    sum_time = 0;
                                    sum_it = 0;
                                end
                            else
                                if mod(i1+i2,15) == 0
                                    print_estimate_transition_prob(total_iterations,index,sum_time/sum_it);
                                    sum_time = 0;
                                    sum_it = 0;
                                end
                            end
                            
                        end
                    end
                end
            end
            
        case '3D'
            
            Nx1 = size(temp_state_partition.X1,2); % Number of points in the first dimension
            Nx2 = size(temp_state_partition.X1,1); % Number of points in the second dimension
            Nx3 = size(temp_state_partition.X1,3); % Number of points in the third dimension
            
            number_of_inputs = length(input_partition); % Number of possible inputs
            
            number_of_MC_simulations = param.MC; % Number of Monte Carlo simulation for each transition
            
            number_of_points = Nx1*Nx2*Nx3;
            temp_values = cell(Nx1*Nx2*Nx3*number_of_inputs,1);
            
            total_iterations = Nx1*Nx2*Nx3*number_of_inputs; % total number of remaining iterations
            PrintEstimateTransitionProb(total_iterations,0,0);
            
            % Variable to estimate the amount of remaining time
            sum_time = 0;
            sum_it = 0;
            
            
            for i1 = 1:Nx1
                for i2 = 1:Nx2
                    for i3 = 1:Nx3
                        x = [temp_state_partition.X1(i2,i1,i3);temp_state_partition.X2(i2,i1,i3);temp_state_partition.X3(i2,i1,i3)]; % Get the current state
                        if number_of_MC_simulations <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                            
                            % Iterate over all possible input combination
                            tempTime = tic;
                            temp_prob = zeros(number_of_points,number_of_inputs);
                            parfor j =1:number_of_inputs
                                u = input_partition(j,:); % selecting a particular allowable input
                                prob_xu = RunMonteCarlo(x,u,obj_state_partition,type_vector_field,param); % empirical estimate of the transition probability
                                temp_prob(:,j) = prob_xu;
                            end
                            
                            index = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],number_of_inputs,[]);
                            
                            for j=number_of_inputs-1:-1:0
                                temp_values{index - j} = temp_prob(:,number_of_inputs-j);
                            end
                            
                            final_time_input = toc(tempTime);
                            sum_time = sum_time + final_time_input;
                            sum_it = sum_it + 1;
                            
                            % The two lines below updates the progress bar
                            index = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],number_of_inputs,[]);
                            
                            % Printing results on the screen
                            if (number_of_points-1)/number_of_inputs < 100
                                if mod(i1+i2+i3,5) == 0
                                    PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                    sum_time = 0;
                                    sum_it = 0;
                                end
                            else
                                if mod(i1+i2+i3,15) == 0
                                    PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                    sum_time = 0;
                                    sum_it = 0;
                                end
                            end
                            
                        else % if number of MC simulation is larger than 100, we will leverage parallel computation
                            
                            % Iterate over all possible input combination
                            sum_time = 0;
                            sum_it = 0;
                            for j =1:number_of_inputs
                                u = input_partition(j,:); % selecting a particular allowable input
                                
                                tempTime1 = tic;
                                
                                prob_xu = RunMonteCarloParallel(x,u,obj_state_partition,type_vector_field,12,param); % empirical estimate of the transition probability using parallel computation
                                
                                index = RemainingIterations(4,[[i1;i2;i3;j],[Nx1;Nx2;Nx3;number_of_inputs]],1,[]);
                                temp_values{index} = prob_xu;
                                
                                final_time_input = toc(tempTime1);
                                sum_time = sum_time + final_time_input;
                                sum_it = sum_it + 1;
                                
                                % Printing on the screen
                                if (number_of_points-1)/number_of_inputs < 100
                                    if mod(i1+i2+i3,5) == 0
                                        PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                        sum_time = 0;
                                        sum_it = 0;
                                    end
                                else
                                    if mod(i1+i2+i3,15) == 0
                                        PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                        sum_time = 0;
                                        sum_it = 0;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
            
        case '4D'
            
            Nx1 = size(temp_state_partition.X1,1); % Number of points in the first dimension
            Nx2 = size(temp_state_partition.X1,2); % Number of points in the second dimension
            Nx3 = size(temp_state_partition.X1,3); % Number of points in the third dimension
            Nx4 = size(temp_state_partition.X1,4); % Number of points in the third dimension
            
            number_of_inputs = length(input_partition); % Number of possible inputs
            
            number_of_MC_simulations = param.MC; % Number of Monte Carlo simulation for each transition
            
            number_of_points = Nx1*Nx2*Nx3*Nx4;
            temp_values = cell(number_of_points*number_of_inputs,1);
            
            total_iterations = number_of_points*number_of_inputs; % total number of remaining iterations
            PrintEstimateTransitionProb(total_iterations,0,0);
            
            % Variable to estimate the amount of remaining time
            sum_time = 0;
            sum_it = 0;
            
            
            for i1 = 1:Nx1
                for i2 = 1:Nx2
                    for i3 = 1:Nx3
                        for i4 =1:Nx4
                            
                            x = [temp_state_partition.X1(i1,i2,i3,i4);temp_state_partition.X2(i1,i2,i3,i4);temp_state_partition.X3(i1,i2,i3,i4);temp_state_partition.X4(i1,i2,i3,i4)]; % Get the current state
                            if number_of_MC_simulations <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                                
                                % Iterate over all possible input combination
                                tempTime = tic;
                                temp_prob = zeros(number_of_points,number_of_inputs);
                                for j =1:number_of_inputs
                                    u = input_partition(j,:); % selecting a particular allowable input
                                    prob_xu = RunMonteCarlo(x,u,obj_state_partition,type_vector_field,param); % empirical estimate of the transition probability
                                    temp_prob(:,j) = prob_xu;
                                end
                                
                                index = RemainingIterations(4,[[i1;i2;i3;i4],[Nx1;Nx2;Nx3;Nx4]],number_of_inputs,[]);
                                
                                for j=number_of_inputs-1:-1:0
                                    temp_values{index - j} = temp_prob(:,number_of_inputs-j);
                                end
                                
                                final_time_input = toc(tempTime);
                                sum_time = sum_time + final_time_input;
                                sum_it = sum_it + 1;
                                
                                % Printing results on the screen
                                if (number_of_points-1)/number_of_inputs < 100
                                    if mod(i1+i2+i3,5) == 0
                                        PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                        sum_time = 0;
                                        sum_it = 0;
                                    end
                                else
                                    if mod(i1+i2+i3,15) == 0
                                        PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                        sum_time = 0;
                                        sum_it = 0;
                                    end
                                end
                                
                            else % if number of MC simulation is larger than 100, we will leverage parallel computation
                                
                                % Iterate over all possible input combination
                                sum_time = 0;
                                sum_it = 0;
                                for j =1:number_of_inputs
                                    u = input_partition(j,:); % selecting a particular allowable input
                                    
                                    tempTime1 = tic;
                                    
                                    prob_xu = RunMonteCarloParallel(x,u,obj_state_partition,type_vector_field,12,param); % empirical estimate of the transition probability using parallel computation
                                    
                                    index = RemainingIterations(5,[[i1;i2;i3;i4;j],[Nx1;Nx2;Nx3;Nx4;number_of_inputs]],1,[]);
                                    temp_values{index} = prob_xu;
                                    
                                    final_time_input = toc(tempTime1);
                                    sum_time = sum_time + final_time_input;
                                    sum_it = sum_it + 1;
                                    
                                    % Printing on the screen
                                    if (number_of_points-1)/number_of_inputs < 100
                                        if mod(i1+i2+i3+i4,5) == 0
                                            PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                            sum_time = 0;
                                            sum_it = 0;
                                        end
                                    else
                                        if mod(i1+i2+i3+i4,15) == 0
                                            PrintEstimateTransitionProb(total_iterations,index,sum_time/sum_it);
                                            sum_time = 0;
                                            sum_it = 0;
                                        end
                                    end
                                    
                                end
                            end
                        end
                    end
                end
            end
            
            
        otherwise
            NotImplemented();
    end
    
    out = containers.Map(grid_with_inputs,temp_values);

end

function out = run_monte_carlo(x,u,state_partition,type_vector_field,param)

% This function computes an estimate of the transition P(.|x,u) associated
% with state partition.
%
%
%   Input: x,u -- state-action pair from which we would like to estimate
%                 the transition
%          state_partition -- is a state partition object 
%          type_vector_field -- Type of the vector field.
%          param -- structure with the following parameters:
%                   number_of_MC_simulations -- number of Monte Carlo simulation
%                   grid_no_inputs -- List with the name of the states of the MDP
%                   all fields required by generateNoise,
%                   computeElementPartition functions defined above

number_of_MC_simulations = param.number_of_MC_simulations; % Number of Monte Carlo simulation
grid_no_inputs = param.grid_no_inputs; % Name of the discrete states of the MDP

temp_partition = state_partition.get_values.partition;


number_of_points = length(grid_no_inputs);

temp_values = repmat({0},1,number_of_points);


switch type_vector_field
    case 'TCL'
        
        N = size(temp_partition.X,1);
        
        for i = 1:number_of_MC_simulations % iterate over the number of simulations
            noise = generate_noise(param,type_vector_field); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = state_partition.compute_element_partition(x,u,noise,param);
            index_member_partition = temp.element_partition;
            
            index = remaining_iterations(1,[index_member_partition,N],1,[]);
            temp_values{index} = temp_values{index} + 1/number_of_MC_simulations;
        end

    case 'LTI'
         Nx1 = size(temp_partition.X1,2); % Number of points in the first dimension
         Nx2 = size(temp_partition.X1,1); % Number of points in the second dimension


         for i = 1:number_of_MC_simulations % iterate over the number of simulations
             noise = generate_noise_LTI(param.mean_noise,param.chol_cov); % simulates a noise transition

             % Compute the associated element of the partition and store the index
             temp = state_partition.compute_element_partition(x,u,noise,param);
             index_member_partition = temp.element_partition;

             index = remaining_iterations(2,[index_member_partition,[Nx1;Nx2]],1,[]); % getting the index of the current state
             temp_values{index} = temp_values{index} + 1/number_of_MC_simulations;
         end
        
%     case 'ChainInt'
%         
%         Nx1 = size(temp_partition.X1,2); % Number of points in the first dimension
%         Nx2 = size(temp_partition.X1,1); % Number of points in the second dimension
%         
%         
%         for i = 1:number_of_MC_simulations % iterate over the number of simulations
%             noise = generateNoise(param,type_vector_field); % simulates a noise transition
%             
%             % Compute the associated element of the partition and store the index
%             temp = state_partition.computeElementPartition(x,u,noise,param);
%             index_member_partition = temp.elementPartition;
%             
%             index = RemainingIterations(2,[index_member_partition,[Nx1;Nx2]],1,[]); % getting the index of the current state
%             temp_values{index} = temp_values{index} + 1/number_of_MC_simulations;
%         end
        
%     case 'Fishery'
%         
%         Nx1 = size(temp_partition.X1,2); % Number of points in the first dimension
%         Nx2 = size(temp_partition.X1,1); % Number of points in the second dimension
%         Nx3 = size(temp_partition.X1,3); % Number of points in the third dimension
%         
%         
%         for i = 1:number_of_MC_simulations % iterate over the number of simulations
%             noise = generateNoise(param,type_vector_field); % simulates a noise transition
%             
%             % Compute the associated element of the partition and store the index
%             temp = state_partition.computeElementPartition(x,u,noise,param);
%             index_member_partition = temp.elementPartition;
%             
%             index = RemainingIterations(3,[index_member_partition,[Nx1;Nx2;Nx3]],1,[]); % getting the index of the current state
%             temp_values{index} = temp_values{index} + 1/number_of_MC_simulations;
%         end
%         
%     case 'CarPole'
%         
%         Nx1 = size(temp_partition.X1,1); % Number of points in the first dimension
%         Nx2 = size(temp_partition.X1,2); % Number of points in the second dimension
%         Nx3 = size(temp_partition.X1,3); % Number of points in the third dimension
%         Nx4 = size(temp_partition.X1,4); % Number of points in the third dimension
%         
%         
%         for i = 1:number_of_MC_simulations % iterate over the number of simulations
%             noise = generateNoise(param,type_vector_field); % simulates a noise transition
%             
%             % Compute the associated element of the partition and store the index
%             temp = state_partition.computeElementPartition(x,u,noise,param);
%             index_member_partition = temp.elementPartition;
%             
%             index = RemainingIterations(4,[index_member_partition,[Nx1;Nx2;Nx3;Nx4]],1,[]); % getting the index of the current state
%             temp_values{index} = temp_values{index} + 1/number_of_MC_simulations;
%         end
%         
%     case 'CarPoleNL'
%         
%         Nx1 = size(temp_partition.X1,1); % Number of points in the first dimension
%         Nx2 = size(temp_partition.X1,2); % Number of points in the second dimension
%         Nx3 = size(temp_partition.X1,3); % Number of points in the third dimension
%         Nx4 = size(temp_partition.X1,4); % Number of points in the third dimension
%         
%         
%         for i = 1:number_of_MC_simulations % iterate over the number of simulations
%             noise = generateNoise(param,type_vector_field); % simulates a noise transition
%             
%             % Compute the associated element of the partition and store the index
%             temp = state_partition.computeElementPartition(x,u,noise,param);
%             index_member_partition = temp.elementPartition;
%             
%             index = RemainingIterations(4,[index_member_partition,[Nx1;Nx2;Nx3;Nx4]],1,[]); % getting the index of the current state
%             temp_values{index} = temp_values{index} + 1/number_of_MC_simulations;
%         end
     
    otherwise
        not_implemented();
        
end

out = zeros(number_of_points,1);

for i=1:length(out)
    out(i) = temp_values{i};
end


end

function out = run_monte_carlo_parallel(current_state,current_input,...
                                            state_partition,type_vector_field,...
                                                n_core,param)

% Exploits parallel computation using the RunMonteCarlo function defined
% above. All the input parameters are defined in RunMonteCarlo definition.

temp_state_partition = state_partition.get_values.partition;
number_of_MC_simulations = param.number_of_MC_simulations; % number of Monte Carlo simulation
temp_MC_parallel = round(number_of_MC_simulations/n_core); % rounding of parallel computation using the number of available cores
MC_parallel = zeros(n_core,1); 

% Defining the number of simulation for each core
for i =1:n_core-1
    MC_parallel(i) = temp_MC_parallel;
end

MC_parallel(end) = number_of_MC_simulations - (n_core-1)*temp_MC_parallel; % the remaining simulations are assigned to the last core

if sum(MC_parallel) ~= number_of_MC_simulations
    error('Error in the partition of the MC simulations'); % outputs an error if an inconsistency is found
end


switch type_vector_field
    
    case 'TCL'
        number_of_points = size(temp_state_partition.X,1);
        
        temp = zeros(number_of_points,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        parfor i = 1:n_core
            temp_param = param;
            temp_param.number_of_MC_simulations = MC_parallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = run_monte_carlo(current_state,current_input,state_partition,type_vector_field,temp_param);
        end
        
        out = zeros(number_of_points,1);

    case 'LTI'
        % Size of the state partition
        Nx1 = size(temp_state_partition.X1,2);
        Nx2 = size(temp_state_partition.X1,1);
        
        temp = zeros(Nx1*Nx2,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        parfor i = 1:n_core
            temp_param = param;
            temp_param.MC = MC_parallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = run_monte_carlo(current_state,current_input,state_partition,type_vector_field,temp_param);
        end
        
        out = zeros(Nx1*Nx2,1);
        
%     case 'ChainInt'
%         
%         % Size of the state partition
%         Nx1 = size(temp_state_partition.X1,2);
%         Nx2 = size(temp_state_partition.X1,1);
%         
%         temp = zeros(Nx1*Nx2,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
%         
%         for i = 1:n_core
%             temp_param = param;
%             temp_param.MC = MC_parallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
%             temp(:,i) = RunMonteCarlo(current_state,current_input,state_partition,type_vector_field,temp_param);
%         end
%         
%         out = zeros(Nx1*Nx2,1);
%     
%     case 'Fishery'
%         % Size of the state partition
%         Nx1 = size(temp_state_partition.X1,2);
%         Nx2 = size(temp_state_partition.X1,1);
%         Nx3 = size(temp_state_partition.X1,3);
%         
%         temp = zeros(Nx1*Nx2*Nx3,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
%         
%         for i = 1:n_core
%             temp_param = param;
%             temp_param.MC = MC_parallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
%             temp(:,i) = RunMonteCarlo(current_state,current_input,state_partition,type_vector_field,temp_param);
%         end
%         
%         out = zeros(Nx1*Nx2*Nx3,1);
%         
%     case 'CarPole'
%         
%         % Size of the state partition
%         Nx1 = size(temp_state_partition.X1,1);
%         Nx2 = size(temp_state_partition.X1,2);
%         Nx3 = size(temp_state_partition.X1,3);
%         Nx4 = size(temp_state_partition.X1,4);
%         
%         temp = zeros(Nx1*Nx2*Nx3*Nx4,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
%         
%         for i = 1:n_core
%             temp_param = param;
%             temp_param.MC = MC_parallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
%             temp(:,i) = RunMonteCarlo(current_state,current_input,state_partition,type_vector_field,temp_param);
%         end
%         
%         out = zeros(Nx1*Nx2*Nx3*Nx4,1);
%         
%     case 'CarPoleNL'
%         
%          % Size of the state partition
%         Nx1 = size(temp_state_partition.X1,1);
%         Nx2 = size(temp_state_partition.X1,2);
%         Nx3 = size(temp_state_partition.X1,3);
%         Nx4 = size(temp_state_partition.X1,4);
%         
%         temp = zeros(Nx1*Nx2*Nx3*Nx4,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
%         
%         for i = 1:n_core
%             temp_param = param;
%             temp_param.MC = MC_parallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
%             temp(:,i) = RunMonteCarlo(current_state,current_input,state_partition,type_vector_field,temp_param);
%         end
%         
%         out = zeros(Nx1*Nx2*Nx3*Nx4,1);
        
        
    otherwise
        NotImplemented();
end


% merging the result of the parallel computation
for i =1:n_core
    out = out + temp(:,i)*(MC_parallel(i)/number_of_MC_simulations);
end


end



