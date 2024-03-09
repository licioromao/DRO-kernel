function [] = add_function_paths(path_project)


addpath(strcat(path_project,'/AmbiguityClasses/'));
addpath(strcat(path_project,'/AuxiliaryFunctions/'));
addpath(strcat(path_project,'/EstimatingTransProb/'));
addpath(strcat(path_project,'/Kernels/'));
addpath(strcat(path_project,'/PartitionClassesAndFunctions/'));
addpath(strcat(path_project,'/PlotFunctions/'));
addpath(strcat(path_project,'/ValueFunctionComputation/'));
addpath(strcat(path_project,'/ValueFunctionComputation/TCL_VectorField'));
addpath(strcat(path_project,'/ValueFunctionComputation/LTI_VectorField'));
addpath(strcat(path_project,'/VectorFields/'));
addpath(genpath(strcat(path_project,'/Dependencies/')));


end