function out = CarPole(TimeHorizon,NumberOfPartitions,NumberOfMonteCarlo,StructAmbiguityTypes,OuterLoopInfo,FILE)

TransitionProb = [];

if nargin > 5
    load(FILE,'TransitionProb')
    load(FILE,'-regexp','ValueFunc*');
else
    FILE = [];
end


%% Model parameters 

T = 0.2; % Discretization 
A = [[0,1,0,0];[0,-10.95,-2.75,0.0043];[0,0,0,1];[0,24.92,28.48,-0.044]];
B = [0;1.94;0;-4.44];

Al = [-0.7;-1;-pi/6;-10*pi]; % Lower bound on the safe set
Ah = -Al; % Upper bound on the safe set

mu = 0; sigma = 0.01^2; % Estimate on the mean and variance

param.T = T;
param.A = A;
param.B = B;
param.SafeSet = [Al,Ah];
param.ReachSet = [[-30;-30;-0.05;-30],[30;30;0.05;30]];

param.N = int8(TimeHorizon);
param.OuterLoopInfo = OuterLoopInfo;
param.OuterLoopInfo.StringAmbiguity = StructAmbiguityTypes; % THIS IS BAD PRACTICE SINCE I AM ALSO SHARING AMBIGUITY PARAMETERS

pd = makedist('Normal','mu',mu,'sigma',sigma);

param.w = pd;

param.MC = NumberOfMonteCarlo;

param.NumberOfPartitions = NumberOfPartitions;
InputPartition = generateInputPartition([],'CarPole'); % Generate a vector with all possible combinations of inputs


Grid = StatePartition(NumberOfPartitions,param.SafeSet,'CarPole'); % Generate the partition of the state space
[List,ListX] = Grid.createList(InputPartition); % List containing labels for the discrete states of the discretazation

param.List = List;
param.ListX = ListX;
param.sizePartition = Grid.getSizePartition;

% Generating the transition probability
if isempty(TransitionProb) % only executes if TransitionPorb is empty
    TransitionProb = EstimateTransition(Grid,InputPartition,param); 
end

param.TransitionProb = TransitionProb;
% Value function computation
% 
% % Saving the transitions
% if nargin > 4
%     save(FILE);
% else
%     FileName = getDateSaveFile(NumberOfPartitions,'TCL'); % Getting the name of file based on the current date and time
%     save(FileName); % saving the results in ./results/
% end



L = length(StructAmbiguityTypes);

for i = 1:L
    switch StructAmbiguityTypes{i}.Name
        case 'NoAmbiguity'
            TotalTime = tic;
            ValueFuncNoAmbiguity = MainValueFunctionIteration(Grid,InputPartition,'CarPole',StructAmbiguityTypes{i},exist('ValueFuncNoAmbiguity','var'),param);
            ValueFuncNoAmbiguity.time = toc(TotalTime);
        case 'MomentAmbiguity'
            TotalTime1 = tic;
            ValueFuncMoment = MainValueFunctionIteration(Grid,InputPartition,'CarPole',StructAmbiguityTypes{i},exist('ValueFuncMoment','var'),param);
            ValueFuncMoment.time = toc(TotalTime1);
            
            paramSave.rhoMu = StructAmbiguityTypes{i}.rhoMu;
            paramSave.rhoSigma = StructAmbiguityTypes{i}.rhoSigma;
             
        case 'WassersteinAmbiguity'
            TotalTime2 = tic;
            ValueFuncWasserstein = MainValueFunctionIteration(Grid,InputPartition,'CarPole',StructAmbiguityTypes{i},exist('ValueFuncWasserstein','var'),param);
            ValueFuncWasserstein.time = toc(TotalTime2);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
                   
            
        case 'KLdivAmbiguity'
            TotalTime3 = tic;
            ValueFuncKL = MainValueFunctionIteration(Grid,InputPartition,'CarPole',StructAmbiguityTypes{i},exist('ValueFuncKL','var'),param);
            ValueFuncKL.time = toc(TotalTime3);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
            
        case 'KernelAmbiguity'  
            TotalTime4 = tic;
            ValueFuncKernel = MainValueFunctionIteration(Grid,InputPartition,'CarPole',StructAmbiguityTypes{i},exist('ValueFuncKernel','var'),param);
            ValueFuncKernel.time = toc(TotalTime4);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
           
        otherwise
            warning('%s has not been implemented. Jumping to the next string',StructAmbiguityTypes{i});
    end   
end

if nargin > 5
    save(FILE);
else
    FileName = getDateSaveFile(TimeHorizon,NumberOfPartitions,NumberOfMonteCarlo,'CarPole',paramSave); % Getting the name of file based on the current date and time
    save(FileName.FullPath); % saving the results in the path specified by FILE
end

out = FileName;


end