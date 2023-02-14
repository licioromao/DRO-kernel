function [] = remove_function_paths(path_project)

rmpath(strcat(path_project,'/AmbiguityClasses/'));
%rmpath(strcat(path_project,'/AuxiliaryFunctions/'));
rmpath(strcat(path_project,'/EstimatingTransProb/'));
rmpath(strcat(path_project,'/Kernels/'));
rmpath(strcat(path_project,'/PartitionClassesAndFunctions/'));
rmpath(strcat(path_project,'/PlotFunctions/'));
rmpath(strcat(path_project,'/ValueFunctionComputation/'));
rmpath(strcat(path_project,'/VectorFields/'));

rmpath(genpath(strcat(path_project,'/Dependencies/')));

end