% 12/08/2022, licio.romao@gmail.com
% 
% This code implements the example in Section 6.1 of paper "Verification of
% discrete time stochastic hybrid systems: A stochastic reach-avoid
% decision problem" by Summers and Lygeros, Automatica, 2010.

% The code uses the functions defined at the bottom of this file as well as
% the BLABLA



function out = Fishery_management(TimeHorizon,NumberOfPartitions,NumberOfMonteCarlo,StructAmbiguityTypes,FILE)

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
InputPartition = generateInputPartition(Input,'Fishery'); % Generate a vector with all possible combinations of inputs
param.SizeInputPartition = size(InputPartition,1);

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

Grid = StatePartition(NumberOfPartitions,param.SafeSet,'Fishery'); % Generate the partition of the state space
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

L = length(StructAmbiguityTypes);

for i = 1:L
    switch StructAmbiguityTypes{i}.Name
        case 'NoAmbiguity'
            TotalTime = tic;
            ValueFuncNoAmbiguity = MainValueFunctionIteration(Grid,InputPartition,'Fishery',StructAmbiguityTypes{i},exist('ValueFuncNoAmbiguity','var'),param);
            ValueFuncNoAmbiguity.time = toc(TotalTime);
        case 'MomentAmbiguity'
            TotalTime1 = tic;
            ValueFuncMoment = MainValueFunctionIteration(Grid,InputPartition,'Fishery',StructAmbiguityTypes{i},exist('ValueFuncMoment','var'),param);
            ValueFuncMoment.time = toc(TotalTime1);
            
            paramSave.rhoMu = StructAmbiguityTypes{i}.rhoMu;
            paramSave.rhoSigma = StructAmbiguityTypes{i}.rhoSigma;
        case 'WassersteinAmbiguity'
            TotalTime2 = tic;
            ValueFuncWasserstein = MainValueFunctionIteration(Grid,InputPartition,'Fishery',StructAmbiguityTypes{i},exist('ValueFuncWasserstein','var'),param);
            ValueFuncWasserstein.time = toc(TotalTime2);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
            
        case 'KLdivAmbiguity'
            TotalTime3 = tic;
            ValueFuncKL = MainValueFunctionIteration(Grid,InputPartition,'Fishery',StructAmbiguityTypes{i},exist('ValueFuncKL','var'),param);
            ValueFuncKL.time = toc(TotalTime3);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
            
        case 'KernelAmbiguity'  
            TotalTime4 = tic;
            ValueFuncKernel = MainValueFunctionIteration(Grid,InputPartition,'Fishery',StructAmbiguityTypes{i},exist('ValueFuncKernel','var'),param);
            ValueFuncKernel.time = toc(TotalTime4);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
            
        otherwise
            warning('%s has not been implemented. Jumping to the next string',StructAmbiguityTypes{i});
    end
end

if nargin > 4
    save(FILE);
else
    FileName = getDateSaveFile(TimeHorizon,NumberOfPartitions,NumberOfMonteCarlo,'Fishery',paramSave); % Getting the name of file based on the current date and time
    save(FileName.FullPath); % saving the results in ./results/
end

out = FileName;

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
