function out = TCL(TimeHorizon,NumberOfPartitions,NumberOfMonteCarlo,AmbiguityTypes,FILE)

TransitionProb = [];

if nargin > 4
    load(FILE,'TransitionProb')
    load(FILE,'-regexp','ValueFunc*');
else
    FILE = [];
end


%% Model parameters 

R = 2; % Thermal resistence in Celsius/kW
C = 2; % Thermal capacity ini kW/Celsius
P = 14; % Range of energy transfer to or from the thermal mass in kW
eta = 0.7; % Control efficiency
h = 5/60; % Discretization time in hours

alpha = exp(-h/(C*R));

Al = 19; % Lower bound on the safe set
Ah = 22; % Upper bound on the safe set

mu = 0; sigma = 0.25^2; % Estimate on the mean and variance
W = 50*[-0.5*sqrt(sigma/12),0.5*sqrt(sigma/12)];  % Support of the distribution

theta = 32; % Environment temperature
p = 0.05; % success probability

param.alpha = alpha;
param.theta = theta;
param.eta = eta;
param.R = R;
param.P = P;
param.C = C;
param.h = h;
param.p = p;
param.SafeSet = [Al,Ah];
param.N = int8(TimeHorizon);

pd = makedist('Normal','mu',mu,'sigma',sigma);
pd_truncated = truncate(pd,W(1),W(2));

param.w = pd_truncated;

param.MC = NumberOfMonteCarlo;

param.NumberOfPartitions = NumberOfPartitions;
InputPartition = generateInputPartition([],'TCL'); % Generate a vector with all possible combinations of inputs


Grid = StatePartition(NumberOfPartitions,param.SafeSet,'TCL'); % Generate the partition of the state space
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
            ValueFuncNoAmbiguity = MainValueFunctionIteration(Grid,InputPartition,'TCL',AmbiguityTypes{i},exist('ValueFuncNoAmbiguity','var'),param);
        case 'MomentAmbiguity'
            ValueFuncMoment = MainValueFunctionIteration(Grid,InputPartition,'TCL',AmbiguityTypes{i},exist('ValueFuncMoment','var'),param);
        case 'WassersteinAmbiguity'
            ValueFuncWasserstein = MainValueFunctionIteration(Grid,InputPartition,'TCL',AmbiguityTypes{i},exist('ValueFuncWasserstein','var'),param);
        case 'KLdivAmbiguity'
            ValueFuncKL = MainValueFunctionIteration(Grid,InputPartition,'TCL',AmbiguityTypes{i},exist('ValueFuncKL','var'),param);
        case 'KernelAmbiguity'  
            ValueFuncKernel = MainValueFunctionIteration(Grid,InputPartition,'TCL',AmbiguityTypes{i},exist('ValueFuncKernel','var'),param);
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

