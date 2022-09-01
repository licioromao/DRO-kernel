% 12/08/2022, licio.romao@gmail.com
% 
% This code implements the example in Section 6.1 of paper "Verification of
% discrete time stochastic hybrid systems: A stochastic reach-avoid
% decision problem" by Summers and Lygeros, Automatica, 2010.

% The code uses the functions defined at the bottom of this file as well as
% the BLABLA

%% Model parameters 

N = 10; % Time horizon of the reach-avoid property
L = 200; % Biomass limit
r = 1; % per-capita recruitment 
C = 70; % Maximum target catch

No_MonteCarlo = 50; % Number of BLABLA
    
% All distributions are considered Gaussian with the parameters given as
% below
mu_v = 0.2; sigma_v = 0.1^2;
mu_gamma = 1; sigma_gamma = 0.6^2;
mu_lambda = 0.25; sigma_lambda = 0.1^2;
mu_delta = 1.1; sigma_delta = 0.2^2;

v = makedist('Normal','mu',mu_v,'sigma',sigma_v);
lambda = makedist('Normal','mu',mu_lambda,'sigma',sigma_lambda);
delta = makedist('Normal','mu',mu_delta,'sigma',sigma_delta);
gamma = makedist('Normal','mu',mu_gamma,'sigma',sigma_gamma);


Input = [0;0.5;1]; % Discretized control space. In this case, the formulation of the problem only allows for three distinct control inputs
InputPartition = generateInputPartition(Input); % Generate a vector with all possible combinations of inputs
NumberInputs = size(InputPartition,1); % Number of all possible combination of control inputs

SafeSet = [0, 200;
           0, 200;
           0, 1200]; % Safe set
       
ReachSet = [100, 200; 
            100, 200; 
            800, 1200]; % Reach set

% All the parameters of the model will be stored in the following structure
        
param.N = N;
param.L = L;
param.r = r;
param.C = C;
param.SafeSet = SafeSet;
param.ReachSet = ReachSet;

param.delta = delta;
param.v = v;
param.lambda = lambda;
param.gamma = gamma;
param.MC = No_MonteCarlo;

%%% Dynamical simulation

% The implementation of the dynamics uses three user-defined functions that
% can be found in the end of this .m file
 
% x = zeros(N_trajectories,N_steps + 1); 
% input = zeros(N_trajectories,N_steps);
% noise = zeros(N_trajectories,N_steps);
% 
% for ell = 1:N_trajectories
%      temp = generateTrajec(@f,x0,N_steps,param);
%      x(ell,:) = temp.x;
%      input(ell,:) = temp.u;
%      noise(ell,:) = temp.noise;
% end


% plot(X_scale,x) % Plot the generated trajectories

%% Partitioning the state space

No_partition = [3,3,3]; % Define the number of points at each dimensional of the state space
param.NumberOfPartitions = No_partition;
%Noise = generateNoise(param);

Grid = StatePartition(No_partition,param.SafeSet); % Generate the partition of the state space
List = Grid.createList; % List containing labels for the discrete states of the discretazation

param.List = List;
param.sizePartition = Grid.getSizePartition;

%% Generating the transition probability
TransitionProb = EstimateTransition(Grid,InputPartition,param); 

param.TransitionProb = TransitionProb;

%% Value function computation

NumberOfPoints = No_partition(1)*No_partition(2)*No_partition(3); % number of points of the value function
valueFunction = zeros(NumberOfPoints + 1,N+1); % initializing the variable to store the value function
OptInput = zeros(NumberOfPoints + 1,N+1); % initializing the variable to store the optimal input

IndexSafeAndReachSet = getIndexReachAvoid(Grid.getValues,param); % getting index of the safe and reach set

param.IndexSafeAndReachSet = IndexSafeAndReachSet; % adding the index to the param structure

valueFunction(IndexSafeAndReachSet.reachIndex,N+1) = 1; % initializing the value function on the reach set

final = 1;

% Creating a progress bar of the value function computation
h = waitbar(0,'Initializing','Name','Computing Value Function...');
total_iterations = (N-final)*NumberOfPoints*NumberInputs;
fprintf('Total of iterations: %d \n',total_iterations);

for i = N:-1:final
    NextValueFunc = valueFunction(:,i+1); % saving in a temporary variable the value function of the next step
    parfor j = 1:NumberOfPoints % iterates over the number of points
        x = Grid.getValues.grid_x(j,:)'; % getting the current state to be updated
        tempValueFunc = zeros(NumberInputs,1);
        for uCounter = 1:NumberInputs
            u = InputPartition(uCounter,:)'; % iterating over the number of inputs
            tempValueFunc(uCounter) = iterateValueFunction(NextValueFunc,x,u,param); % getting the new value for the value function
        end
        [valueFunction(j,i),OptInput(j,i)] = max(tempValueFunc); % storing the optimal value function and policy
    end
    
    % updating the progress bar
    index = (N-i)*NumberOfPoints*NumberInputs + NumberOfPoints*NumberInputs;
    waitbar(index/total_iterations,h,sprintf('%.5f completed',index/total_iterations));
end

%delete(h)
%save
%% Below we try to plot the results for given values of the trid variable

close all

InitialState_X3 = 750;
InitialState_X3 = Grid.getValues.X3(1,1,find(Grid.getValues.X3(1,1,:) <= InitialState_X3,1,'last'));

[XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid.getValues,param);

levels = [0.1,0.2,0.4,0.5,0.6,0.8,0.85,0.9,1];

figure 
h = contour(XX,YY,ProjValueFunc,levels,'ShowText','on','LineWidth',2);


figure
pcolor(XX,YY,OptPolicy)


% 
% InitialState_X3 = 420;
% InitialState_X3 = Grid.X3(1,1,find(Grid.X3(1,1,:) <= InitialState_X3,1,'last'));
% 
% [XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid,param);
% 
% levels = [0.1,0.2,0.4,0.5,0.6,0.8,0.85,0.9,1];
% 
% figure
% h = contour(XX,YY,ProjValueFunc,levels,'ShowText','on','LineWidth',2);
% 
% 
% figure
% pcolor(XX,YY,OptPolicy)




% 
% InitialState_X3 = 600;
% InitialState_X3 = Grid.X3(1,1,find(Grid.X3(1,1,:) <= InitialState_X3,1,'last'));
% 
% [XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid,param);
% 
% figure 
% h = contour(XX,YY,ProjValueFunc,levels,'ShowText','on','LineWidth',2);

%% Functions related to the partition and generating an estimate of the transition kernel

function out = generateInputPartition(Input)
% This function generates a grid of size Nu^2 x 2 that contains any
% possible combination of discrete control inputs. The number of column is
% equal to the number of patches, which in this example is equal to 2.
    
    Nu = length(Input);
    out = zeros(Nu^2,2); % Initializing the output
    
    for i=1:Nu
        for j =1:Nu
            out(j+Nu*(i-1),1) = Input(i);
            out(j+Nu*(i-1),2) = Input(j);
        end
    end

end

%% Computing the transition probabilities

function out = RunMonteCarlo(x,u,Grid,param)

% This function computes an estimate of the transition P(.|x,u) associated
% with state partition.
%
%
%   Input: func -- vector field function handle. See definition of the
%                  function above
%          x,u -- state-action pair from which we would like to estimate
%                 the transition
%          StatePartition -- represents the output of the function
%                            generatePartition above 
%          param -- structure with the following parameters:
%                   MC -- number of Monte Carlo simulation
%                   List -- List with the name of the states of the MDP
%                   all fields required by generateNoise,
%                   computeElementPartition functions defined above
%    This function also uses the Hashtable class function

StatePartition = Grid.getValues;

MC = param.MC; % Number of Monte Carlo simulation
List = param.List; % Name of the discrete states of the MDP
hTable = HashTable([]); % Create a hashtable. See the corresponding file for additional information

hTable.List = List; % Initialize the hashtable with the name of the states of the MDP
hTable.lengthList = length(List); % stores the number of states

for i = 1:MC % iterate over the number of simulations
    noise = generateNoise(param); % simulates a noise transition
    
    % Compute the associated element of the partition and store the index
    temp = Grid.computeElementPartition(x,u,noise,param); 
    indexMemberPartition = temp.elementPartition;
    
    if ~isempty(indexMemberPartition) % if an element if found
        
        memberPartition = [StatePartition.X1(1,indexMemberPartition(1),1),StatePartition.X2(indexMemberPartition(2),1,1),StatePartition.X3(1,1,indexMemberPartition(3))]; % get the corresponding element
        
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

StatePartition = Grid.getValues;
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
Nx1 = size(StatePartition.X1,2); 
Nx2 = size(StatePartition.X1,1);
Nx3 = size(StatePartition.X1,3);

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

function out = EstimateTransition(Grid,InputPartition,param)

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

StatePartition = Grid.getValues;

Nx1 = size(StatePartition.X1,2); % Number of points in the first dimension
Nx2 = size(StatePartition.X1,1); % Number of points in the second dimension
Nx3 = size(StatePartition.X1,3); % Number of points in the third dimension

Nu = length(InputPartition); % Number of possible inputs

MC = param.MC; % Number of Monte Carlo simulation for each transition

out = []; 

h = waitbar(0,'Initializing','Name','Generating Transition probabilities...'); % Plotting a bar to check the progress
total_iterations = Nx1*Nx2*Nx3*Nu; % total number of remaining iterations
n_core = 4; % number of cores. AJUST THIS PARAMETER ACCORDING TO YOUR COMPUTATIONAL POWER

for i1 = 1:Nx1
    for i2 = 1:Nx2
        for i3 = 1:Nx3
            tic; % time estimate for one iteration
            x = [StatePartition.X1(i2,i1,i3),StatePartition.X2(i2,i1,i3),StatePartition.X3(i2,i1,i3)]; % Get the current state
            stringX = sprintf('(%.2f,%.2f,%.2f)',x(1),x(2),x(3)); % Convert the state into an allowable name
            temp = cell(Nu,1); % temp variable
            if MC <= 100 % if number of MC simulations if less than 100 do not use parallel computation
                
                % Iterate over all possible input combination
                parfor j =1:Nu
                    u = InputPartition(j,:); % selecting a particular allowable input
                    stringU = sprintf('(%.2f,%.2f)',u(1),u(2)); % creating the corresponding string
                    
                    prob_xu = RunMonteCarlo(x,u,Grid,param); % empirical estimate of the transition probability
                    
                    % saving the results
                    temp{j}.State = stringX; 
                    temp{j}.Action = stringU;
                    temp{j}.ProbMeasure = prob_xu;
                end
                
                % The two lines below updates the progress bar
                index = (i1-1)*Nx2*Nx3*Nu + (i2-1)*Nx3*Nu + (i3-1)*Nu + Nu; 
                waitbar(index/total_iterations,h,sprintf('%.5f completed',index/total_iterations));
                
            else % if number of MC simulation is larger than 100, we will leverage parallel computation
                
                 % Iterate over all possible input combination
                for j =1:Nu
                    u = InputPartition(j,:); % selecting a particular allowable input
                    stringU = sprintf('(%.2f,%.2f)',u(1),u(2)); % creating the corresponding string
                    
                    prob_xu = RunMonteCarloParallel(x,u,Grid,100,param); % empirical estimate of the transition probability using parallel computation
                    
                    % saving results
                    temp{j}.State = stringX;
                    temp{j}.Action = stringU;
                    temp{j}.ProbMeasure = prob_xu;
                    
                    % The two lines below updates the progress bar
                    index = (i1-1)*Nx2*Nx3*Nu + (i2-1)*Nx3*Nu + (i3-1)*Nu + j;
                    waitbar(index/total_iterations,h,sprintf('%.5f completed',index/total_iterations));
                end
            end
            estimate_time = toc; % storing time for a single state-action pair estimate 
            (estimate_time*(Nx1*Nx2*Nx3 - (i1-1)*Nx2*Nx3 - (i2-1)*Nx3 - i3))/(60*60*24) % estimate of the remaining time in days
            out = [out;temp];
        end
    end
end

delete(h) % deleting the progress bar

end

%% Computation of the value function

function out = getIndexReachAvoid(StatePartition,param)

% Retuns the indices of the safe and reach sets associated with the partition

Nx1 = size(StatePartition.X1,2); 
Nx2 = size(StatePartition.X1,1);
Nx3 = size(StatePartition.X1,3);

% Information about the safe and reach sets
SafeSet = param.SafeSet;
ReachSet = param.ReachSet;

% Initializing the variables that will store the indices
safeIndex = false(Nx1*Nx2*Nx3 +1,1);
reachIndex = false(Nx1*Nx2*Nx3 + 1,1);

% Iterating over states
for i1 = 1:Nx1
    for i2=1:Nx2
        for i3=1:Nx3
            x = [StatePartition.X1(i2,i1,i3);StatePartition.X2(i2,i1,i3);StatePartition.X3(i2,i1,i3)];  % current state
            index = (i1-1)*Nx2*Nx3 + (i2-1)*Nx3 + i3;
            
            % if current state belongs to safe set, then set the
            % corresponding index to true
            if min(SafeSet(:,1) <= x) && min(x <= SafeSet(:,2)) 
                safeIndex(index) = true;
            else
                safeIndex(index) = false;
            end
            
            % do the same with the reach set
            if min(ReachSet(:,1) <= x) && min(x <= ReachSet(:,2))
                reachIndex(index) = true;
            else
                reachIndex(index) = false;
            end
        end
    end
end

% The last state correspond to the unsafe state of the MDP, so we set to
% false
safeIndex(end) = false;
reachIndex(end) = false;


% saving the results
out.safeIndex = safeIndex;
out.reachIndex = reachIndex;

end

function out = valueFunctionH(x,u,Z,param)

% Returns the value function for a given state-action pair (x,u) using the
% value function at the next iteration (Z). Param is a structure with
% fields necessary for all the subroutines (check each subroutine for this precise information)
    
TransitionProb = param.TransitionProb; % vector with the transition probability matrix
L = length(TransitionProb); % number of states in the transition probability

if length(Z) ~= length(TransitionProb{1}.ProbMeasure)
    error('The dimension of the thrid argument is incorrect.') % outputs an error if L is inconsistent with the size of the input Z
end

stringX = sprintf('(%.2f,%.2f,%.2f)',x(1),x(2),x(3)); % creating the string of the current state
stringU = sprintf('(%.2f,%.2f)',u(1),u(2)); % creating the string of the current action

key = false;
i = 1;

% Searching for the correct transition probability
while i <= L && ~key
    if strcmp(TransitionProb{i}.State,stringX) && strcmp(TransitionProb{i}.Action,stringU) % if stringX and stringU matches with the information in the transition prob matrix
        key = true;
        out = TransitionProb{i}.ProbMeasure'*Z; % value function at the current state-action pair
    end
    i = i+1;
end

if ~key % outputs an error is there is no element in TransitionProb with the label (stringX,stringU)
    error('The input-action pair is not a member of the transition probability')
end

end

function out = iterateValueFunction(ValueFunc,CurrentState,Input,param)

indexCurrentState = find(strcmp(param.List,sprintf('(%.2f,%.2f,%.2f)',CurrentState(1),CurrentState(2),CurrentState(3))));

indexReachSet = param.IndexSafeAndReachSet.reachIndex;

if indexReachSet(indexCurrentState)
    out = 1;
else
    out = valueFunctionH(CurrentState,Input,ValueFunc,param);
end

end

%% Functions created to assist the plotting of the obtained results
function [X,Y,ProjValueFunction,OptPolicy] = plot_results(Height,ValueFunction,OptInput,...
            StatePartition,param)

Indices = find(StatePartition.grid_x(:,3) == Height);

if isempty(Indices)
    error('The required initial condition is not compatible with the state partition');
else
    gridXY = StatePartition.grid_x(:,1:2);
    x = unique(gridXY(:,1));
    y = unique(gridXY(:,2));
    [X,Y] = meshgrid(x,y);
    
    tempValue = ValueFunction(Indices);
    tempPolicy = OptInput(Indices);
    
    ProjValueFunction = reshape(tempValue,size(X));
    OptPolicy = reshape(tempPolicy,size(X));
end


end
















