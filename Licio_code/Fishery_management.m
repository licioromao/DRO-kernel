% 12/08/2022, licio.romao@gmail.com
% 
% This code implements the example in Section 6.1 of paper "Verification of
% discrete time stochastic hybrid systems: A stochastic reach-avoid
% decision problem" by Summers and Lygeros, Automatica, 2010.

% The code uses the functions defined at the bottom of this file as well as
% the BLABLA



function out = Fishery_management(TimeHorizon,NumberOfPartitions,NumberOfMonteCarlo,AmbiguityTypes,FILE)

TransitionProb = [];

if nargin > 4
    load(FILE,'TransitionProb')
    load(FILE,'-regexp','ValueFunc*');
else
    FILE = [];
end

% Model parameters 
L = 200; % Biomass limit
r = 1; % per-capita recruitment
C = 70; % Maximum target catch

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

param.N = TimeHorizon;
param.L = L;
param.r = r;
param.C = C;
param.SafeSet = SafeSet;
param.ReachSet = ReachSet;

param.delta = delta;
param.v = v;
param.lambda = lambda;
param.gamma = gamma;
param.MC = NumberOfMonteCarlo;

% Partitioning the state space

param.NumberOfPartitions = NumberOfPartitions;
%Noise = generateNoise(param);

Grid = StatePartition(NumberOfPartitions,param.SafeSet); % Generate the partition of the state space
List = Grid.createList; % List containing labels for the discrete states of the discretazation

param.List = List;
param.sizePartition = Grid.getSizePartition;

% Generating the transition probability
if isempty(TransitionProb) % only executes if TransitionPorb is empty
    TransitionProb = EstimateTransition(Grid,InputPartition,param);
end

param.TransitionProb = TransitionProb;
% Value function computation

L = length(AmbiguityTypes);

for i = 1:L
    switch AmbiguityTypes{i}
        case 'NoAmbiguity'
            ValueFuncNoAmbiguity = MainValueFunctionIteration(Grid,InputPartition,AmbiguityTypes{i},exist('ValueFuncNoAmbiguity','var'),param);
        case 'MomentAmbiguity'
            ValueFuncMoment = MainValueFunctionIteration(Grid,InputPartition,AmbiguityTypes{i},exist('ValueFuncMoment','var'),param);
        case 'WassersteinAmbiguity'
            ValueFuncWasserstein = MainValueFunctionIteration(Grid,InputPartition,AmbiguityTypes{i},exist('ValueFuncWasserstein','var'),param);
        case 'KLdivAmbiguity'
            ValueFuncKL = MainValueFunctionIteration(Grid,InputPartition,AmbiguityTypes{i},exist('ValueFuncKL','var'),param);
        case 'KernelAmbiguity'  
            ValueFuncKernel = MainValueFunctionIteration(Grid,InputPartition,AmbiguityTypes{i},exist('ValueFuncKernel','var'),param);
        otherwise
            warning('%s has not been implemented. Jumping to the next string',AmbiguityTypes{i});
    end
end

if nargin > 4
    save(FILE);
else
    FileName = getDateSaveFile(NumberOfPartitions); % Getting the name of file based on the current date and time
    save(FileName); % saving the results in ./results/
end

out = [];

end

% Below we try to plot the results for given values of the trid variable

% close all
% 
% InitialState_X3 = 750;
% InitialState_X3 = Grid.getValues.X3(1,1,find(Grid.getValues.X3(1,1,:) <= InitialState_X3,1,'last'));
% 
% [XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid.getValues,param);
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

%Functions related to the partition and generating an estimate of the transition kernel

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

% Computing the transition probabilities

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
                    index = RemainingIterations(3,[[i1;i2;i3;j],[Nx1;Nx2;Nx3;Nu]],1,h);
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

% Functions created to assist the plotting of the obtained results
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


% Function that computes the value function

function ValueFunc = MainValueFunctionIteration(PartitionObj,InputPartition,TypeOfAmbiguity,BooleanSignal,param)

if ~BooleanSignal
    switch TypeOfAmbiguity
        case 'NoAmbiguity'
            fprintf('\nComputing value function (without ambiguity)...\n')
        case 'MomentAmbiguity'
            fprintf('\nComputing value function (moment ambiguity)...\n')
        case 'KernelAmbiguity'
            fprintf('\nComputing value function (kernel ambiguity)...\n')
        case 'KLdivAmbiguity'
            fprintf('\nComputing value function (KL ambiguity)...\n')
        case 'WassersteinAmbiguity'
            fprintf('\nComputing value function (Wasserstein ambiguity)...\n')
        otherwise
            error('Type of Ambiguity not implemented')
    end
    
    tic;
    ValueFunc = ComputeValueFunction(param,TypeOfAmbiguity);
    
    ValueFunc = ValueFunc.getIndexReachAvoid(PartitionObj.getValues);
    ValueFunc = ValueFunc.BackwardIteration(PartitionObj,InputPartition);
    ValueFunc.time = toc;
    
    fprintf('Done\n');
else
    ValueFunc = [];
    switch TypeOfAmbiguity
        case 'NoAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncNoAmbiguity');
        case 'MomentAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncMoment');
        case 'KernelAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncKernel');
        case 'KLdivAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncKLdivAmbiguity');
        case 'WassersteinAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncWasserstein');
        otherwise
            error('Type of Ambiguity not implemented')
    end
end








end













