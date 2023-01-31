function [] = AddFunctionPaths(pathProject)


addpath(strcat(pathProject,'/AmbiguityClasses/'));
addpath(strcat(pathProject,'/AuxiliaryFunctions/'));
addpath(strcat(pathProject,'/EstimatingTransProb/'));
addpath(strcat(pathProject,'/Kernels/'));
addpath(strcat(pathProject,'/PartitionClassesAndFunctions/'));
addpath(strcat(pathProject,'/PlotFunctions/'));
addpath(strcat(pathProject,'/ValueFunctionComputation/'));
addpath(strcat(pathProject,'/VectorFields/'));
addpath(genpath(strcat(pathProject,'/Dependencies/')));


end