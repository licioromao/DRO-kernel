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

No_MonteCarlo = 800; % Number of BLABLA
    
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

No_partition = [20,20,33]; % Define the number of points at each dimensional of the state space
%Noise = generateNoise(param);

Grid = generatePartition(No_partition,param.SafeSet); % Generate the partition of the state space
List = createList(Grid); % List containing labels for the discrete states of the discretazation

param.List = List;
param.sizePartition = getSizePartition(Grid);

%% Generating the transition probability
TransitionProb = EstimateTransition(@vector_field,Grid,InputPartition,param); 

param.TransitionProb = TransitionProb;

%% Value function computation

NumberOfPoints = No_partition(1)*No_partition(2)*No_partition(3); % number of points of the value function
valueFunction = zeros(NumberOfPoints + 1,N+1); % initializing the variable to store the value function
OptInput = zeros(NumberOfPoints + 1,N+1); % initializing the variable to store the optimal input

IndexSafeAndReachSet = getIndexReachAvoid(Grid,param); % getting index of the safe and reach set

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
        x = Grid.grid_x(j,:)'; % getting the current state to be updated
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
InitialState_X3 = Grid.X3(1,1,find(Grid.X3(1,1,:) <= InitialState_X3,1,'last'));

[XX,YY,ProjValueFunc,OptPolicy] = plot_results(InitialState_X3,valueFunction(:,1),OptInput(:,1),Grid,param);

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


%% User-defined functions that define the vector field

function out = vector_field_R(x,param)

% This function implements the function R that composes the vector field,
% as given in the first column of page 1956 of the paper mentioned at the
% top of this file.

% Inputs:   x - current state
%           param - this is a structure that contains the field r and L

r = param.r;
L = param.L;

out = max(0,r*x*(1-(x/L))); % output of the function
end

function [out,mode] = vector_field_C(x,d,delta,param)

% This function implements the function C as shown in page 1956 of the
% paper at the top of this file.

%  Inputs: x - current state
%          d - represents the input of the dynamics
%          delta - is the variability of the catch
%          param -  is a structure containing L and C as fields
%
%  Output: out - value of the expression
%          mode - mode of the system. Note that this system has two modes,
%          according to the value of the current state and L
%

L = param.L;
C = param.C;

if x < L % if x is less than L
    out = delta*d*C*x/L;
    mode = 1;
else % otherwise
    out = delta*d*C;
    mode = 2;
end

end

function out = vector_field_E(x1,x2,d,v,delta,param)

% This function implements the function E as shown in page 1956 of the
% paper at the top of this file.

% Inputs: x1,x2 -  fish biomass in each patch
%         d - represents the input of the dynamics
%         delta - is the variability of the catch
%         v - natural fish mortality
%         param - is a structure containing L and C as fields

out = (1-v)*x1 - vector_field_C(x1+x2,d,delta,param); % output of this expression as shown in the paper

end

function out = vector_field_M(x1,x2,d1,d2,v1,v2,delta1,delta2,param)

% This function implements the function M as shown in page 1956 of the
% paper at the top of this file, which represents the fish escapament
% between patch at a given period

% Inputs: x1,x2 -  fish biomass in each patch
%         d1,d2 - represents the inputs for each patch 
%         delta1,delta2 - is the catch variability in each patch
%         v1,v2 - natural fish mortality in each patch
%         param - is a structure containing L and C as fields


temp_1 = vector_field_E(x1,x2,d1,v1,delta1,param);
temp_2 = vector_field_E(x2,x1,d2,v2,delta2,param);

out = temp_1 - temp_2;
end

function out = generateNoise(param)

% This function generates the noise that acts on the dynamics of the system.

% Input: param - structure containing the distributions associated with any
% source of noise in the systems

out = [random(param.v,2,1);random(param.gamma,2,1);random(param.lambda);random(param.delta,2,1)];

end

function out = vector_field(CurrentState,Input,Noise,param)

% Equation (14)-(16) of the paper mentioned at the top of this file.

% Input: Noise - a 7-dimensional vector containing the sources of noise in
%        the model
%        Input - Control input   
%        CurrentState - current state of the model
%        param - is a structure that contains all the parameters of the
%        model
%
% Output: out - next state of the model


v1 = Noise(1);
v2 = Noise(2);
gamma1 = Noise(3);
gamma2 = Noise(4);
lambda = Noise(5);
delta1 = Noise(6);
delta2 = Noise(7);

R1 = vector_field_R(CurrentState(1),param);
R2 = vector_field_R(CurrentState(2),param);

M = vector_field_M(CurrentState(1),CurrentState(2),Input(1),Input(2),v1,v2,delta1,delta2,param);

[C1,mode] = vector_field_C(CurrentState(1)+CurrentState(2),Input(1),delta1,param);
[C2,~] = vector_field_C(CurrentState(1)+CurrentState(2),Input(2),delta2,param);

out = [(1-v1)*CurrentState(1) + gamma1*R1 - C1 - lambda*M; % Equation (14)
        (1-v2)*CurrentState(2) + gamma2*R2 - C2 + lambda*M; % Equation (15)
        CurrentState(3) + C1 + C2; % Equation (16)
        mode]; % Mode of the system

end

%% Functions related to the partition and generating an estimate of the transition kernel

function out = generatePartition(No_partition,SafeSet)

% Generates a lattice partition of the state space, taking as limit the
% safe set
%
% Input: No_Partition - a vector containing the number of points at each
%                       state of the model
%
%        SafeSet - A 3x2 matrix with the limits of the safe set in each
%                  dimension
%
% Output: out - is a structure with four fields as defined below:
%               out.X1, out.X2, out.X3 -- multi-dimensional matrices of size No_Partition(1)
%                                         x No_Partition(2) x No_Partition(3) containing the points
%                                         of the generated lattice
%               out.grid_x -- elements of the grid organized in a (No_Partition(1)
%                                         x No_Partition(2) x
%                                         No_Partition(3)) x 3 matrix
% The field out.grid_x is used to assign labels to the discrete states of
% the model.

grid_x1 = linspace(SafeSet(1,1),SafeSet(1,2),No_partition(1));
grid_x2 = linspace(SafeSet(2,1),SafeSet(2,2),No_partition(2));
grid_x3 = linspace(SafeSet(3,1),SafeSet(3,2),No_partition(3));

[X1,X2,X3] = meshgrid(grid_x1,grid_x2,grid_x3);

Nx1 = size(X1,2);
Nx2 = size(X1,1);
Nx3 = size(X1,3);

grid_x = zeros(Nx1*Nx2*Nx3,3);

for i1 = 1:Nx1
    for i2 = 1:Nx2
        for i3 = 1:Nx3
            x = [X1(i2,i1,i3);X2(i2,i1,i3);X3(i2,i1,i3)];
            index = (i1-1)*Nx2*Nx3 + (i2-1)*Nx3 + i3;
            grid_x(index,:) = x;
        end
    end
end

out.X1 = X1;
out.X2 = X2;
out.X3 = X3;
out.grid_x = grid_x;

end

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

function out = createList(StatePartition)

% This function creates a vector of array that contains the labels of the
% states of the generated discrete model

% Inputs: StatePartition -- this is the output of the function
% generatePartition as defined above
%
% Ouput: out -- labels for the discrete states

Nx1 = size(StatePartition.X1,2);
Nx2 = size(StatePartition.X2,1);
Nx3 = size(StatePartition.X3,3);

out = {};

for i1 = 1:Nx1
    for i2 = 1:Nx2
        for i3 = 1:Nx3
            x = [StatePartition.X1(i2,i1,i3),StatePartition.X2(i2,i1,i3),StatePartition.X3(i2,i1,i3)];
            out = [out;sprintf('(%.2f,%.2f,%.2f)',x(1),x(2),x(3))];
        end
    end
end

out = [out;'NaN']; % The last state called 'NaN' represents a fictious state of the discrete model containing unsafe states

end

function out = getSizePartition(Partition)

% This function returns the granularity of the partition in each dimension.
% Its input is the output of the function generatePartition.

    h1 = Partition.X1(1,2,1) - Partition.X1(1,1,1);
    h2 = Partition.X2(2,1,1) - Partition.X2(1,1,1);
    h3 = Partition.X3(1,1,2) - Partition.X3(1,1,1);
    
    out = [h1;h2;h3];
end

function out = getElementPartition(x,Partition)

% This function returns the element of the partition that x belongs to. If
% the input x is outside the safe set, this function returns an empty
% array.

 index = [find(Partition.X1(1,:,1) <= x(1),1,'last');
        find(Partition.X2(:,1,1) <= x(2),1,'last');
        find(Partition.X3(1,1,:) <= x(3),1,'last')]; % Check the element of the partition
    
 if length(index) ~= 3 % If there is one dimensional for which index is empty, set out to empty array
     index = [];
     out = [];
 end
 
 if ~isempty(index) % Otherwise, the function outputs a structure with index and xHat field
     xHat = [Partition.X1(1,index(1),1);Partition.X2(index(2),1,1);Partition.X3(1,1,index(3));x(4)];
     out.index = index;
     out.x = xHat;
 end
 
end

function out = getCenterPartition(x,sizePartition)

% This function receives as input a point x in the lattice and returns the
% center of the partition.

out = zeros(4,1);

out(1:3) = x(1:3) + sizePartition/2;
out(4) = x(4);

end

function out = computeElementPartition(func,currentState,input,noise,Partition,param)

% Given currentState, input and noise vectors, this function uses the
% vector_field function to progate the dynamics and returns the element and
% center of the partition corresponding to the next state.

out.nextState = [];
out.elementPartition = [];

if ~isempty(func)
    nextState = func(currentState,input,noise,param);
    temp = getElementPartition(nextState,Partition);
    out.elementPartition = temp.index;
    out.nextState = getCenterPartition(temp.x,param.sizePartition);
else
    out.elementPartition = getElementPartition(currentState,Partition);
end

end

%% Computing the transition probabilities

function out = RunMonteCarlo(func,x,u,StatePartition,param)

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


MC = param.MC; % Number of Monte Carlo simulation
List = param.List; % Name of the discrete states of the MDP
hTable = HashTable([]); % Create a hashtable. See the corresponding file for additional information

hTable.List = List; % Initialize the hashtable with the name of the states of the MDP
hTable.lengthList = length(List); % stores the number of states

for i = 1:MC % iterate over the number of simulations
    noise = generateNoise(param); % simulates a noise transition
    
    % Compute the associated element of the partition and store the index
    temp = computeElementPartition(func,x,u,noise,StatePartition,param); 
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

function out = RunMonteCarloParallel(func,x,u,StatePartition,n_core,param)

% Exploits parallel computation using the RunMonteCarlo function defined
% above. All the input parameters are defined in RunMonteCarlo definition.

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
    temp(:,i) = RunMonteCarlo(func,x,u,StatePartition,tempParam);
end

out = zeros(Nx1*Nx2*Nx3 + 1,1);

% merging the result of the parallel computation
for i =1:n_core
    out = out + temp(:,i)*(MCParallel(i)/MC);
end

end

function out = EstimateTransition(func,StatePartition,InputPartition,param)

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
                    
                    prob_xu = RunMonteCarlo(func,x,u,StatePartition,param); % empirical estimate of the transition probability
                    
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
                    
                    prob_xu = RunMonteCarloParallel(func,x,u,StatePartition,100,param); % empirical estimate of the transition probability using parallel computation
                    
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
















