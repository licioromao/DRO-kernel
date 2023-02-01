function [] = remove_function_paths(pathProject)

rmpath(strcat(pathProject,'/AmbiguityClasses/'));
rmpath(strcat(pathProject,'/AuxiliaryFunctions/'));
rmpath(strcat(pathProject,'/EstimatingTransProb/'));
rmpath(strcat(pathProject,'/Kernels/'));
rmpath(strcat(pathProject,'/PartitionClassesAndFunctions/'));
rmpath(strcat(pathProject,'/PlotFunctions/'));
rmpath(strcat(pathProject,'/ValueFunctionComputation/'));
rmpath(strcat(pathProject,'/VectorFields/'));

rmpath(genpath(strcat(pathProject,'/Dependencies/')));

end