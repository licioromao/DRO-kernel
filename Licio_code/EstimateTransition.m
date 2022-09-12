
function out = EstimateTransition(Grid,InputPartition,TypeOfVectorFields,param)

% Iterates over all state-action pair and compute the transition
% probability.
%
%
%   Input: func -- handle of the vector field that defines the dynamics
%          StatePartition -- output of the generatePartition function
%          InputPartition -- output of the generateInput partition function
%          param -- a structure with all the field required by the
%          subfunctions that are needed (please refer to the function definition above or to the main code for an example of an allowable structure)
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

StatePartition = Grid.getValues.Partition;

Nx1 = size(StatePartition.X1,2); % Number of points in the first dimension
Nx2 = size(StatePartition.X1,1); % Number of points in the second dimension
Nx3 = size(StatePartition.X1,3); % Number of points in the third dimension

Nu = length(InputPartition); % Number of possible inputs

MC = param.MC; % Number of Monte Carlo simulation for each transition

out = []; 

h = waitbar(0,'Initializing','Name','Generating Transition probabilities...'); % Plotting a bar to check the progress
total_iterations = Nx1*Nx2*Nx3*Nu; % total number of remaining iterations

for i1 = 1:Nx1
    for i2 = 1:Nx2
        for i3 = 1:Nx3
            tic; % time estimate for one iteration
            x = [StatePartition.X1(i2,i1,i3),StatePartition.X2(i2,i1,i3),StatePartition.X3(i2,i1,i3)]; % Get the current state
            stringX = sprintf('(%.2f,%.2f,%.2f)',x(1),x(2),x(3)); % Convert the state into an allowable name
            temp = cell(Nu,1); % temp variable
            if MC <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                
                % Iterate over all possible input combination
                Lastime = 0;
                tic;
                parfor j =1:Nu
                    u = InputPartition(j,:); % selecting a particular allowable input
                    stringU = sprintf('(%.2f,%.2f)',u(1),u(2)); % creating the corresponding string
                    
                    prob_xu = RunMonteCarlo(x,u,Grid,param); % empirical estimate of the transition probability
                    
                    % saving the results
                    temp{j}.State = stringX; 
                    temp{j}.Action = stringU;
                    temp{j}.ProbMeasure = prob_xu;
                end
                Lastime = toc;
                
                % The two lines below updates the progress bar
                index = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],Nu,h); 
                SecToGo = (total_iterations - index)*Lastime;
                waitbar(index/total_iterations,h,sprintf('%.5f completed. %.2f seconds to go.',index/total_iterations,SecToGo));
                
            else % if number of MC simulation is larger than 100, we will leverage parallel computation
                
                 % Iterate over all possible input combination
                for j =1:Nu
                    u = InputPartition(j,:); % selecting a particular allowable input
                    stringU = sprintf('(%.2f,%.2f)',u(1),u(2)); % creating the corresponding string
                    
                    Lastime = 0;
                    tic;
                    prob_xu = RunMonteCarloParallel(x,u,Grid,100,param); % empirical estimate of the transition probability using parallel computation
                    Lastime = toc;
                    % saving results
                    temp{j}.State = stringX;
                    temp{j}.Action = stringU;
                    temp{j}.ProbMeasure = prob_xu;
                    
                    % The two lines below updates the progress bar
                    index = RemainingIterations(4,[[i1;i2;i3;j],[Nx1;Nx2;Nx3;Nu]],1,h);
                    SecToGo = (total_iterations - index)*Lastime;
                    waitbar(index/total_iterations,h,sprintf('%.5f completed. %.2f seconds to go.',index/total_iterations,SecToGo));
                end
            end
            %estimate_time = toc; % storing time for a single state-action pair estimate 
            %(estimate_time*(Nx1*Nx2*Nx3 - (i1-1)*Nx2*Nx3 - (i2-1)*Nx3 - i3))/(60*60*24) % estimate of the remaining time in days
            out = [out;temp];
        end
    end
end

delete(h) % deleting the progress bar

end


% Computing the transition probabilities

function out = RunMonteCarlo(x,u,Grid,param)

% This function computes an estimate of the transition P(.|x,u) associated
% with state partition.
%
%
%   Input: x,u -- state-action pair from which we would like to estimate
%                 the transition
%          Grid -- is a State partition object 
%          param -- structure with the following parameters:
%                   MC -- number of Monte Carlo simulation
%                   List -- List with the name of the states of the MDP
%                   all fields required by generateNoise,
%                   computeElementPartition functions defined above

TypeOfVectorField = Grid.getValues.TypeOfVectorField; % Type of vector field
MC = param.MC; % Number of Monte Carlo simulation
List = param.List; % Name of the discrete states of the MDP
hTable = HashTable([]); % Create a hashtable. See the corresponding file for additional information

hTable.List = List; % Initialize the hashtable with the name of the states of the MDP
hTable.lengthList = length(List); % stores the number of states

tempPartition = Grid.getValues.Partition;

for i = 1:MC % iterate over the number of simulations
    noise = generateNoise(param,TypeOfVectorField); % simulates a noise transition
    
    % Compute the associated element of the partition and store the index
    temp = Grid.computeElementPartition(x,u,noise,param); 
    indexMemberPartition = temp.elementPartition;
    
    if ~isempty(indexMemberPartition) % if an element if found
        
        memberPartition = [tempPartition.X1(1,indexMemberPartition(1),1),tempPartition.X2(indexMemberPartition(2),1,1),tempPartition.X3(1,1,indexMemberPartition(3))]; % get the corresponding element
        
        stringPartition = sprintf('(%.2f,%.2f,%.2f)',memberPartition(1),memberPartition(2),memberPartition(3)); % string to be added to the private list in the hash table
        
        hTable = hTable.appendString(stringPartition); % append the element to the partition, nothing is done if an element is already a member of the list
        hTable = hTable.addValue(stringPartition,1); % add one to the corresponding value
    else % if the transition happens to a state outside the safe region
        stringPartition = 'NaN'; 
        hTable = hTable.appendString(stringPartition); % append the element to the partition, nothing is done if an element is already a member of the list
        hTable = hTable.addValue(stringPartition,1); % add one to the corresponding value
    end
            
end

out = hTable.createProbMeasure(); % computes the emperical probability distribution over discrete states of the MDP

end

function out = RunMonteCarloParallel(x,u,Grid,n_core,param)

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

% Size of the state partition
Nx1 = size(tempStatePartition.X1,2); 
Nx2 = size(tempStatePartition.X1,1);
Nx3 = size(tempStatePartition.X1,3);

temp = zeros(Nx1*Nx2*Nx3 + 1,n_core); % initilizating the variable that contains the intermediate estimates for the transition probabilities

parfor i = 1:n_core
    tempParam = param;
    tempParam.MC = MCParallel(i); % assigning the correct parameter to be passed onto RunMonteCarlo
    temp(:,i) = RunMonteCarlo(x,u,Grid,tempParam);
end

out = zeros(Nx1*Nx2*Nx3 + 1,1);

% merging the result of the parallel computation
for i =1:n_core
    out = out + temp(:,i)*(MCParallel(i)/MC);
end

end
