function out = EstimateTransition(Grid,InputPartition,param)

% This function iterates over all state-action pairs and computes the
% transition probability.
%
%
%   Input: Grid -- This is a StateParition object.
%          InputPartition -- output of the generateInput partition function
%          TypeOfVectorFields -- This is the type of vector fields.
%          Currently we have implemented TCL and Fishery.
%          param -- a structure with all parameters of the problem
%
%
%   Output: TBD


TypeOfVectorField = Grid.getValues.TypeOfVectorField;

List = param.List;

switch TypeOfVectorField
    
    case 'TCL'
        
        out = EstimateProbTransition(Grid,InputPartition,List,'1D',TypeOfVectorField,param);
        
    case 'ChainInt'
        
        out = EstimateProbTransition(Grid,InputPartition,List,'2D',TypeOfVectorField,param);
    
    case 'Fishery'
        
        out = EstimateProbTransition(Grid,InputPartition,List,'3D',TypeOfVectorField,param);
    
    case 'CarPole'
        
        out = EstimateProbTransition(Grid,InputPartition,List,'4D',TypeOfVectorField,param);
        
    case 'CarPoleNL'
        
        out = EstimateProbTransition(Grid,InputPartition,List,'4D',TypeOfVectorField,param);
        
    otherwise
        NotImplemented();
end




end


% Computing the transition probabilities

function out = EstimateProbTransition(Grid,InputPartition,List,sizeSpace,TypeOfVectorField,param)
    

    tempStatePartition = Grid.getValues.Partition;
    
    switch sizeSpace
        case '1D'
            
            N = size(tempStatePartition.X,1); % Number of points in the first dimension
            
            Nu = length(InputPartition); % Number of possible inputs
            
            MC = param.MC; % Number of Monte Carlo simulation for each transition
            
            NumberOfPoints = N;
            tempValues = cell(N*Nu,1);
            
            total_iterations = N*Nu; % total number of remaining iterations
            PrintEstimateTransitionProb(total_iterations,0,0);
            
            
            % Variables to estimate the ramining time
            sumTime = 0;
            sumIt = 0;
            
            
            for i = 1:N
                x = tempStatePartition.X(i); % Get the current state
                if MC <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                    
                    % Iterate over all possible inputs
                    tempTime2 = tic;
                    tempParForProb = zeros(NumberOfPoints,Nu);
                    
                    for j =1:Nu
                        u = InputPartition(j,:); % selecting a particular allowable input
                        prob_xu = RunMonteCarlo(x,u,Grid,TypeOfVectorField,param); % empirical estimate of the transition probability
                        tempParForProb(:,j) = prob_xu;
                    end
                    
                    indexTrans = RemainingIterations(1,[i,N],Nu,[]);
                    
                    for j=Nu-1:-1:0
                        tempValues{indexTrans - j} = tempParForProb(:,Nu-j);
                    end
                    
                    Lastime = toc(tempTime2);
                    sumTime = sumTime + Lastime;
                    sumIt = sumIt + 1;
                    
                    % The two lines below updates the progress bar
                    index = RemainingIterations(1,[i,N],Nu,[]);
                    
                    % Printing on the screen
                    if (NumberOfPoints-1)/Nu < 100
                        if mod(i+j,5) == 0
                            PrintEstimateTransitionProb(total_iterations,index,sumTime/sumIt);
                            sumTime = 0;
                            sumIt = 0;
                        end
                    else
                        if mod(i+j,15) == 0
                            PrintEstimateTransitionProb(total_iterations,index,sumTime/sumIt);
                            sumTime = 0;
                            sumIt = 0;
                        end
                    end
                    
                else % if number of MC simulation is larger than 100, we will leverage parallel computation
                    
                    % Iterate over all possible input combination
                    for j =1:Nu
                        u = InputPartition(j); % selecting a particular allowable input
                        
                        tempTime3 = tic;
                        
                        prob_xu = RunMonteCarloParallel(x,u,Grid,TypeOfVectorField,12,param); % empirical estimate of the transition probability using parallel computation
                        
                        indexTrans = RemainingIterations(2,[[i;j],[N;Nu]],1,[]);
                        tempValues{indexTrans} = prob_xu;
                        
                        Lastime = toc(tempTime3);
                        sumTime = sumTime + Lastime;
                        sumIt = sumIt + 1;
                        
                        % Printing on the screen
                        if (NumberOfPoints-1)/Nu < 100
                            if mod(i+j,5) == 0
                                PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                sumTime = 0;
                                sumIt = 0;
                            end
                        else
                            if mod(i+j,15) == 0
                                PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                sumTime = 0;
                                sumIt = 0;
                            end
                        end
                    end
                end
            end
            
            
        case '2D'
            
            Nx1 = size(tempStatePartition.X1,2); % Number of points in the first dimension
            Nx2 = size(tempStatePartition.X1,1); % Number of points in the second dimension
            
            Nu = length(InputPartition); % Number of possible inputs
            
            MC = param.MC; % Number of Monte Carlo simulation for each transition
            
            NumberOfPoints = Nx1*Nx2;
            tempValues = cell(Nx1*Nx2*Nu,1);
            
            total_iterations = Nx1*Nx2*Nu; % total number of remaining iterations
            PrintEstimateTransitionProb(total_iterations,0,0);
            
            % Variable to estimate the amount of remaining time
            sumTime = 0;
            sumIt = 0;
            
            
            for i1 = 1:Nx1
                for i2 = 1:Nx2
                    
                    x = [tempStatePartition.X1(i2,i1);tempStatePartition.X2(i2,i1)]; % Get the current state
                    if MC <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                        
                        % Iterate over all possible input combination
                        tempTime = tic;
                        tempParForProb = zeros(NumberOfPoints,Nu);
                        parfor j =1:Nu
                            u = InputPartition(j,:); % selecting a particular allowable input
                            prob_xu = RunMonteCarlo(x,u,Grid,TypeOfVectorField,param); % empirical estimate of the transition probability
                            tempParForProb(:,j) = prob_xu;
                        end
                        
                        indexTrans = RemainingIterations(2,[[i1;i2],[Nx1;Nx2]],Nu,[]);
                        
                        for j=Nu-1:-1:0
                            tempValues{indexTrans - j} = tempParForProb(:,Nu-j);
                        end
                        
                        Lastime = toc(tempTime);
                        sumTime = sumTime + Lastime;
                        sumIt = sumIt + 1;
                        
                        % Printing results on the screen
                        if (NumberOfPoints-1)/Nu < 100
                            if mod(i1+i2,5) == 0
                                PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                sumTime = 0;
                                sumIt = 0;
                            end
                        else
                            if mod(i1+i2,15) == 0
                                PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                sumTime = 0;
                                sumIt = 0;
                            end
                        end
                        
                    else % if number of MC simulation is larger than 100, we will leverage parallel computation
                        
                        % Iterate over all possible input combination
                        sumTime = 0;
                        sumIt = 0;
                        for j =1:Nu
                            u = InputPartition(j,:); % selecting a particular allowable input
                            
                            tempTime1 = tic;
                            
                            prob_xu = RunMonteCarloParallel(x,u,Grid,TypeOfVectorField,12,param); % empirical estimate of the transition probability using parallel computation
                            
                            indexTrans = RemainingIterations(3,[[i1;i2;j],[Nx1;Nx2;Nu]],1,[]);
                            tempValues{indexTrans} = prob_xu;
                            
                            Lastime = toc(tempTime1);
                            sumTime = sumTime + Lastime;
                            sumIt = sumIt + 1;
                            
                            % Printing on the screen
                            if (NumberOfPoints-1)/Nu < 100
                                if mod(i1+i2,5) == 0
                                    PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            else
                                if mod(i1+i2,15) == 0
                                    PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            end
                            
                        end
                    end
                end
            end
            
        case '3D'
            
            Nx1 = size(tempStatePartition.X1,2); % Number of points in the first dimension
            Nx2 = size(tempStatePartition.X1,1); % Number of points in the second dimension
            Nx3 = size(tempStatePartition.X1,3); % Number of points in the third dimension
            
            Nu = length(InputPartition); % Number of possible inputs
            
            MC = param.MC; % Number of Monte Carlo simulation for each transition
            
            NumberOfPoints = Nx1*Nx2*Nx3;
            tempValues = cell(Nx1*Nx2*Nx3*Nu,1);
            
            total_iterations = Nx1*Nx2*Nx3*Nu; % total number of remaining iterations
            PrintEstimateTransitionProb(total_iterations,0,0);
            
            % Variable to estimate the amount of remaining time
            sumTime = 0;
            sumIt = 0;
            
            
            for i1 = 1:Nx1
                for i2 = 1:Nx2
                    for i3 = 1:Nx3
                        x = [tempStatePartition.X1(i2,i1,i3);tempStatePartition.X2(i2,i1,i3);tempStatePartition.X3(i2,i1,i3)]; % Get the current state
                        if MC <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                            
                            % Iterate over all possible input combination
                            tempTime = tic;
                            tempParForProb = zeros(NumberOfPoints,Nu);
                            parfor j =1:Nu
                                u = InputPartition(j,:); % selecting a particular allowable input
                                prob_xu = RunMonteCarlo(x,u,Grid,TypeOfVectorField,param); % empirical estimate of the transition probability
                                tempParForProb(:,j) = prob_xu;
                            end
                            
                            indexTrans = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],Nu,[]);
                            
                            for j=Nu-1:-1:0
                                tempValues{indexTrans - j} = tempParForProb(:,Nu-j);
                            end
                            
                            Lastime = toc(tempTime);
                            sumTime = sumTime + Lastime;
                            sumIt = sumIt + 1;
                            
                            % The two lines below updates the progress bar
                            index = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],Nu,[]);
                            
                            % Printing results on the screen
                            if (NumberOfPoints-1)/Nu < 100
                                if mod(i1+i2+i3,5) == 0
                                    PrintEstimateTransitionProb(total_iterations,index,sumTime/sumIt);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            else
                                if mod(i1+i2+i3,15) == 0
                                    PrintEstimateTransitionProb(total_iterations,index,sumTime/sumIt);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            end
                            
                        else % if number of MC simulation is larger than 100, we will leverage parallel computation
                            
                            % Iterate over all possible input combination
                            sumTime = 0;
                            sumIt = 0;
                            for j =1:Nu
                                u = InputPartition(j,:); % selecting a particular allowable input
                                
                                tempTime1 = tic;
                                
                                prob_xu = RunMonteCarloParallel(x,u,Grid,TypeOfVectorField,12,param); % empirical estimate of the transition probability using parallel computation
                                
                                indexTrans = RemainingIterations(4,[[i1;i2;i3;j],[Nx1;Nx2;Nx3;Nu]],1,[]);
                                tempValues{indexTrans} = prob_xu;
                                
                                Lastime = toc(tempTime1);
                                sumTime = sumTime + Lastime;
                                sumIt = sumIt + 1;
                                
                                % Printing on the screen
                                if (NumberOfPoints-1)/Nu < 100
                                    if mod(i1+i2+i3,5) == 0
                                        PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                        sumTime = 0;
                                        sumIt = 0;
                                    end
                                else
                                    if mod(i1+i2+i3,15) == 0
                                        PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                        sumTime = 0;
                                        sumIt = 0;
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
            
        case '4D'
            
            Nx1 = size(tempStatePartition.X1,1); % Number of points in the first dimension
            Nx2 = size(tempStatePartition.X1,2); % Number of points in the second dimension
            Nx3 = size(tempStatePartition.X1,3); % Number of points in the third dimension
            Nx4 = size(tempStatePartition.X1,4); % Number of points in the third dimension
            
            Nu = length(InputPartition); % Number of possible inputs
            
            MC = param.MC; % Number of Monte Carlo simulation for each transition
            
            NumberOfPoints = Nx1*Nx2*Nx3*Nx4;
            tempValues = cell(NumberOfPoints*Nu,1);
            
            total_iterations = NumberOfPoints*Nu; % total number of remaining iterations
            PrintEstimateTransitionProb(total_iterations,0,0);
            
            % Variable to estimate the amount of remaining time
            sumTime = 0;
            sumIt = 0;
            
            
            for i1 = 1:Nx1
                for i2 = 1:Nx2
                    for i3 = 1:Nx3
                        for i4 =1:Nx4
                            
                            x = [tempStatePartition.X1(i1,i2,i3,i4);tempStatePartition.X2(i1,i2,i3,i4);tempStatePartition.X3(i1,i2,i3,i4);tempStatePartition.X4(i1,i2,i3,i4)]; % Get the current state
                            if MC <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                                
                                % Iterate over all possible input combination
                                tempTime = tic;
                                tempParForProb = zeros(NumberOfPoints,Nu);
                                for j =1:Nu
                                    u = InputPartition(j,:); % selecting a particular allowable input
                                    prob_xu = RunMonteCarlo(x,u,Grid,TypeOfVectorField,param); % empirical estimate of the transition probability
                                    tempParForProb(:,j) = prob_xu;
                                end
                                
                                indexTrans = RemainingIterations(4,[[i1;i2;i3;i4],[Nx1;Nx2;Nx3;Nx4]],Nu,[]);
                                
                                for j=Nu-1:-1:0
                                    tempValues{indexTrans - j} = tempParForProb(:,Nu-j);
                                end
                                
                                Lastime = toc(tempTime);
                                sumTime = sumTime + Lastime;
                                sumIt = sumIt + 1;
                                
                                % Printing results on the screen
                                if (NumberOfPoints-1)/Nu < 100
                                    if mod(i1+i2+i3,5) == 0
                                        PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                        sumTime = 0;
                                        sumIt = 0;
                                    end
                                else
                                    if mod(i1+i2+i3,15) == 0
                                        PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                        sumTime = 0;
                                        sumIt = 0;
                                    end
                                end
                                
                            else % if number of MC simulation is larger than 100, we will leverage parallel computation
                                
                                % Iterate over all possible input combination
                                sumTime = 0;
                                sumIt = 0;
                                for j =1:Nu
                                    u = InputPartition(j,:); % selecting a particular allowable input
                                    
                                    tempTime1 = tic;
                                    
                                    prob_xu = RunMonteCarloParallel(x,u,Grid,TypeOfVectorField,12,param); % empirical estimate of the transition probability using parallel computation
                                    
                                    indexTrans = RemainingIterations(5,[[i1;i2;i3;i4;j],[Nx1;Nx2;Nx3;Nx4;Nu]],1,[]);
                                    tempValues{indexTrans} = prob_xu;
                                    
                                    Lastime = toc(tempTime1);
                                    sumTime = sumTime + Lastime;
                                    sumIt = sumIt + 1;
                                    
                                    % Printing on the screen
                                    if (NumberOfPoints-1)/Nu < 100
                                        if mod(i1+i2+i3+i4,5) == 0
                                            PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                            sumTime = 0;
                                            sumIt = 0;
                                        end
                                    else
                                        if mod(i1+i2+i3+i4,15) == 0
                                            PrintEstimateTransitionProb(total_iterations,indexTrans,sumTime/sumIt);
                                            sumTime = 0;
                                            sumIt = 0;
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
    
    out = containers.Map(List,tempValues);

end

function out = RunMonteCarlo(x,u,Grid,TypeOfVectorField,param)

% This function computes an estimate of the transition P(.|x,u) associated
% with state partition.
%
%
%   Input: x,u -- state-action pair from which we would like to estimate
%                 the transition
%          Grid -- is a State partition object 
%          TypeOfVectorFields -- Type of the vector field.
%          param -- structure with the following parameters:
%                   MC -- number of Monte Carlo simulation
%                   ListX -- List with the name of the states of the MDP
%                   all fields required by generateNoise,
%                   computeElementPartition functions defined above

MC = param.MC; % Number of Monte Carlo simulation
List = param.ListX; % Name of the discrete states of the MDP

tempPartition = Grid.getValues.Partition;


NumberOfPoints = length(List);

tempValues = repmat({0},1,NumberOfPoints);


switch TypeOfVectorField
    case 'TCL'
        
        N = size(tempPartition.X,1);
        
        for i = 1:MC % iterate over the number of simulations
            noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = Grid.computeElementPartition(x,u,noise,param);
            indexMemberPartition = temp.elementPartition;
            
            index = RemainingIterations(1,[indexMemberPartition,N],1,[]);
            tempValues{index} = tempValues{index} + 1/MC;
        end
        
    case 'ChainInt'
        
        Nx1 = size(tempPartition.X1,2); % Number of points in the first dimension
        Nx2 = size(tempPartition.X1,1); % Number of points in the second dimension
        
        
        for i = 1:MC % iterate over the number of simulations
            noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = Grid.computeElementPartition(x,u,noise,param);
            indexMemberPartition = temp.elementPartition;
            
            index = RemainingIterations(2,[indexMemberPartition,[Nx1;Nx2]],1,[]); % getting the index of the current state
            tempValues{index} = tempValues{index} + 1/MC;
        end
        
    case 'Fishery'
        
        Nx1 = size(tempPartition.X1,2); % Number of points in the first dimension
        Nx2 = size(tempPartition.X1,1); % Number of points in the second dimension
        Nx3 = size(tempPartition.X1,3); % Number of points in the third dimension
        
        
        for i = 1:MC % iterate over the number of simulations
            noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = Grid.computeElementPartition(x,u,noise,param);
            indexMemberPartition = temp.elementPartition;
            
            index = RemainingIterations(3,[indexMemberPartition,[Nx1;Nx2;Nx3]],1,[]); % getting the index of the current state
            tempValues{index} = tempValues{index} + 1/MC;
        end
        
    case 'CarPole'
        
        Nx1 = size(tempPartition.X1,1); % Number of points in the first dimension
        Nx2 = size(tempPartition.X1,2); % Number of points in the second dimension
        Nx3 = size(tempPartition.X1,3); % Number of points in the third dimension
        Nx4 = size(tempPartition.X1,4); % Number of points in the third dimension
        
        
        for i = 1:MC % iterate over the number of simulations
            noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = Grid.computeElementPartition(x,u,noise,param);
            indexMemberPartition = temp.elementPartition;
            
            index = RemainingIterations(4,[indexMemberPartition,[Nx1;Nx2;Nx3;Nx4]],1,[]); % getting the index of the current state
            tempValues{index} = tempValues{index} + 1/MC;
        end
        
    case 'CarPoleNL'
        
        Nx1 = size(tempPartition.X1,1); % Number of points in the first dimension
        Nx2 = size(tempPartition.X1,2); % Number of points in the second dimension
        Nx3 = size(tempPartition.X1,3); % Number of points in the third dimension
        Nx4 = size(tempPartition.X1,4); % Number of points in the third dimension
        
        
        for i = 1:MC % iterate over the number of simulations
            noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = Grid.computeElementPartition(x,u,noise,param);
            indexMemberPartition = temp.elementPartition;
            
            index = RemainingIterations(4,[indexMemberPartition,[Nx1;Nx2;Nx3;Nx4]],1,[]); % getting the index of the current state
            tempValues{index} = tempValues{index} + 1/MC;
        end
     
    otherwise
        NotImplemented();
        
end

out = zeros(NumberOfPoints,1);

for i=1:length(out)
    out(i) = tempValues{i};
end


end

function out = RunMonteCarloParallel(x,u,Grid,TypeOfVectorField,n_core,param)

% Exploits parallel computation using the RunMonteCarlo function defined
% above. All the input parameters are defined in RunMonteCarlo definition.

tempStatePartition = Grid.getValues.Partition;
MC = param.MC; % number of Monte Carlo simulation
tempMCParallel = round(MC/n_core); % rounding of parallel computation using the number of available cores
MCParallel = zeros(n_core,1); 

% Defining the number of simulation for each core
for i =1:n_core-1
    MCParallel(i) = tempMCParallel;
end

MCParallel(end) = param.MC - (n_core-1)*tempMCParallel; % the remaining simulations are assigned to the last core

if sum(MCParallel) ~= MC
    error('Error in the partition of the MC simulations'); % outputs an error if an inconsistency is found
end


switch TypeOfVectorField
    
    case 'TCL'
        N = size(tempStatePartition.X,1);
        
        temp = zeros(N,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        parfor i = 1:n_core
            tempParam = param;
            tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = RunMonteCarlo(x,u,Grid,TypeOfVectorField,tempParam);
        end
        
        out = zeros(N,1);
        
    case 'ChainInt'
        
        % Size of the state partition
        Nx1 = size(tempStatePartition.X1,2);
        Nx2 = size(tempStatePartition.X1,1);
        
        temp = zeros(Nx1*Nx2,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        for i = 1:n_core
            tempParam = param;
            tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = RunMonteCarlo(x,u,Grid,TypeOfVectorField,tempParam);
        end
        
        out = zeros(Nx1*Nx2,1);
    
    case 'Fishery'
        % Size of the state partition
        Nx1 = size(tempStatePartition.X1,2);
        Nx2 = size(tempStatePartition.X1,1);
        Nx3 = size(tempStatePartition.X1,3);
        
        temp = zeros(Nx1*Nx2*Nx3,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        for i = 1:n_core
            tempParam = param;
            tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = RunMonteCarlo(x,u,Grid,TypeOfVectorField,tempParam);
        end
        
        out = zeros(Nx1*Nx2*Nx3,1);
        
    case 'CarPole'
        
        % Size of the state partition
        Nx1 = size(tempStatePartition.X1,1);
        Nx2 = size(tempStatePartition.X1,2);
        Nx3 = size(tempStatePartition.X1,3);
        Nx4 = size(tempStatePartition.X1,4);
        
        temp = zeros(Nx1*Nx2*Nx3*Nx4,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        for i = 1:n_core
            tempParam = param;
            tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = RunMonteCarlo(x,u,Grid,TypeOfVectorField,tempParam);
        end
        
        out = zeros(Nx1*Nx2*Nx3*Nx4,1);
        
    case 'CarPoleNL'
        
         % Size of the state partition
        Nx1 = size(tempStatePartition.X1,1);
        Nx2 = size(tempStatePartition.X1,2);
        Nx3 = size(tempStatePartition.X1,3);
        Nx4 = size(tempStatePartition.X1,4);
        
        temp = zeros(Nx1*Nx2*Nx3*Nx4,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        for i = 1:n_core
            tempParam = param;
            tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = RunMonteCarlo(x,u,Grid,TypeOfVectorField,tempParam);
        end
        
        out = zeros(Nx1*Nx2*Nx3*Nx4,1);
        
        
    otherwise
        NotImplemented();
end


% merging the result of the parallel computation
for i =1:n_core
    out = out + temp(:,i)*(MCParallel(i)/MC);
end


end



