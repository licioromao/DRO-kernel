
function out = EstimateTransition(Grid,InputPartition,param)

% Iterates over all state-action pair and compute the transition
% probability.
%
%
%   Input: Grid -- output of the generatePartition function
%          InputPartition -- output of the generateInput partition function
%          TypeOfVectorFields -- This is the type of vector fields.
%          Currently we have implemented TCL and Fishery.
%          param -- a structure with all the field required by the
%
%
%   Output: out -- a cell array of structure of dimension equal the number of elements
%                  in the partition plus one (unsafe state), where each element contains the
%                  following fields:
%
%                  State -- String with the name of the state
%                  Action -- String identifying the action
%                  ProbMeasure -- an estimate of the transition probability
%                  P(.|x,u)

tempStatePartition = Grid.getValues.Partition;
TypeOfVectorField = Grid.getValues.TypeOfVectorField;

List = param.List;



switch TypeOfVectorField
    case 'Fishery'
        
        Nx1 = size(tempStatePartition.X1,2); % Number of points in the first dimension
        Nx2 = size(tempStatePartition.X1,1); % Number of points in the second dimension
        Nx3 = size(tempStatePartition.X1,3); % Number of points in the third dimension
        
        Nu = length(InputPartition); % Number of possible inputs
        
        MC = param.MC; % Number of Monte Carlo simulation for each transition
        
        NumberOfPoints = Nx1*Nx2*Nx3*Nu + 1;
        tempValues = cell(NumberOfPoints,1);
        
        h = waitbar(0,'Initializing','Name','Generating Transition probabilities...'); % Plotting a bar to check the progress
        total_iterations = Nx1*Nx2*Nx3*Nu; % total number of remaining iterations
        
        for i1 = 1:Nx1
            for i2 = 1:Nx2
                for i3 = 1:Nx3
                    x = [tempStatePartition.X1(i2,i1,i3),tempStatePartition.X2(i2,i1,i3),tempStatePartition.X3(i2,i1,i3)]; % Get the current state
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
                        
                        for j=0:Nu-1
                            tempValues{indexTrans +j} = tempParForProb(:,j+1);
                        end
                                                
                        Lastime = toc(tempTime);
                        
                        % The two lines below updates the progress bar
                        index = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],Nu,h);
                        SecToGo = (total_iterations - index)*Lastime;
                        waitbar(index/total_iterations,h,sprintf('%.5f completed. %.2f seconds to go.',index/total_iterations,SecToGo));
                        
                    else % if number of MC simulation is larger than 100, we will leverage parallel computation
                        
                        % Iterate over all possible input combination
                        for j =1:Nu
                            u = InputPartition(j,:); % selecting a particular allowable input
                            
                            tempTime1 = tic;
                            
                            prob_xu = RunMonteCarloParallel(x,u,Grid,TypeOfVectorField,eval(getenv('NUMBER_OF_PROCESSORS')),param); % empirical estimate of the transition probability using parallel computation
                            
                            indexTrans = RemainingIterations(4,[[i1;i2;i3;j],[Nx1;Nx2;Nx3;Nu]],1,[]);
                            tempValues{indexTrans} = prob_xu;
                            
                            Lastime = toc(tempTime1);
                            % The two lines below updates the progress bar
                            SecToGo = (total_iterations - indexTrans)*Lastime;
                            waitbar(indexTrans/total_iterations,h,sprintf('%.5f completed. %.2f seconds to go.',indexTrans/total_iterations,SecToGo));
                        end
                    end
                end
            end
        end
        
    case 'TCL'
        
        N = size(tempStatePartition.X,1); % Number of points in the first dimension
        
        Nu = length(InputPartition); % Number of possible inputs
        
        MC = param.MC; % Number of Monte Carlo simulation for each transition
        
               
        h = waitbar(0,'Initializing','Name','Generating Transition probabilities...'); % Plotting a bar to check the progress
        total_iterations = N*Nu; % total number of remaining iterations
        
        out = cell(N + 1,1);
        
        
        for i = 1:N
            x = tempStatePartition.X(i); % Get the current state
            stringX = sprintf('(%.5f)',x); % Convert the state into an allowable name
            temp = cell(Nu,1); % temp variable
            if MC <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                
                % Iterate over all possible input combination
                tempTime2 = tic;
                for j =1:Nu
                    u = InputPartition(j,:); % selecting a particular allowable input
                    stringU = sprintf('(%.2f)',u); % creating the corresponding string
                    
                    prob_xu = RunMonteCarlo(x,u,Grid,TypeOfVectorField,param); % empirical estimate of the transition probability
                    
                    % saving the results
                    temp{j}.State = stringX;
                    temp{j}.Action = stringU;
                    temp{j}.ProbMeasure = prob_xu;
                end
                Lastime = toc(tempTime2);
                
                % The two lines below updates the progress bar
                index = RemainingIterations(1,[i,N],Nu,h);
                SecToGo = (total_iterations - index)*Lastime;
                waitbar(index/total_iterations,h,sprintf('%.5f completed. %.2f seconds to go.',index/total_iterations,SecToGo));
                
            else % if number of MC simulation is larger than 100, we will leverage parallel computation
                
                % Iterate over all possible input combination
                for j =1:Nu
                    u = InputPartition(j); % selecting a particular allowable input
                    stringU = sprintf('(%.2f)',u); % creating the corresponding string
                    
                    tempTime3 = tic;
                    
                    prob_xu = RunMonteCarloParallel(x,u,Grid,TypeOfVectorField,eval(getenv('NUMBER_OF_PROCESSORS')),param); % empirical estimate of the transition probability using parallel computation
                    
                    % saving results
                    temp{j}.State = stringX;
                    temp{j}.Action = stringU;
                    temp{j}.ProbMeasure = prob_xu;
                    
                    Lastime = toc(tempTime3);
                    
                    % The two lines below updates the progress bar
                    index = RemainingIterations(2,[[i;j],[N;Nu]],1,h);
                    SecToGo = (total_iterations - index)*Lastime;
                    waitbar(index/total_iterations,h,sprintf('%.5f completed. %.2f seconds to go.',index/total_iterations,SecToGo));
                end
            end
            indexTemp = RemainingIterations(1,[i,N],1,h);
            out{indexTemp} = temp;
        end
end

out = containers.Map(List,tempValues);

delete(h) % deleting the progress bar

end


% Computing the transition probabilities

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
%                   List -- List with the name of the states of the MDP
%                   all fields required by generateNoise,
%                   computeElementPartition functions defined above

MC = param.MC; % Number of Monte Carlo simulation
List = param.ListX; % Name of the discrete states of the MDP

tempPartition = Grid.getValues.Partition;


NumberOfPoints = length(List);

tempValues = repmat({0},1,NumberOfPoints);


switch TypeOfVectorField
    case 'Fishery'
        
        Nx1 = size(tempPartition.X1,2); % Number of points in the first dimension
        Nx2 = size(tempPartition.X1,1); % Number of points in the second dimension
        Nx3 = size(tempPartition.X1,3); % Number of points in the third dimension
        
        
        for i = 1:MC % iterate over the number of simulations
            noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = Grid.computeElementPartition(x,u,noise,param);
            indexMemberPartition = temp.elementPartition;
            
            if ~isempty(indexMemberPartition) % if an element if found
                
                index = RemainingIterations(3,[indexMemberPartition,[Nx1;Nx2;Nx3]],1,[]);
                
                tempValues{index} = tempValues{index} + 1/MC;
            else % if the transition happens to a state outside the safe region
                tempValues{end} = tempValues{end} + 1/MC;
            end
            
        end
        
    case 'TCL'
        
        for i = 1:MC % iterate over the number of simulations
            noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
            
            % Compute the associated element of the partition and store the index
            temp = Grid.computeElementPartition(x,u,noise,param);
            indexMemberPartition = temp.elementPartition;
            
            if ~isempty(indexMemberPartition) % if an element if found
                
                memberPartition = tempPartition.X(indexMemberPartition); % get the corresponding element
                
                stringPartition = sprintf('(%.5f)',memberPartition); % string to be added to the private list in the hash table
                
                hTable = hTable.appendString(stringPartition); % append the element to the partition, nothing is done if an element is already a member of the list
                hTable = hTable.addValue(stringPartition,1); % add one to the corresponding value
            else % if the transition happens to a state outside the safe region
                stringPartition = 'NaN';
                hTable = hTable.appendString(stringPartition); % append the element to the partition, nothing is done if an element is already a member of the list
                hTable = hTable.addValue(stringPartition,1); % add one to the corresponding value
            end
            
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
    case 'Fishery'
        % Size of the state partition
        Nx1 = size(tempStatePartition.X1,2);
        Nx2 = size(tempStatePartition.X1,1);
        Nx3 = size(tempStatePartition.X1,3);
        
        temp = zeros(Nx1*Nx2*Nx3 + 1,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        parfor i = 1:n_core
            tempParam = param;
            tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = RunMonteCarlo(x,u,Grid,TypeOfVectorField,tempParam);
        end
        
        out = zeros(Nx1*Nx2*Nx3 + 1,1);
               
    case 'TCL'
        N = size(tempStatePartition.X,1);
        
        temp = zeros(N + 1,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities
        
        for i = 1:n_core
            tempParam = param;
            tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
            temp(:,i) = RunMonteCarlo(x,u,Grid,TypeOfVectorField,tempParam);
        end
        
        out = zeros(N + 1,1);
        
    otherwise
        NotImplemented();
end


% merging the result of the parallel computation
for i =1:n_core
    out = out + temp(:,i)*(MCParallel(i)/MC);
end


end
