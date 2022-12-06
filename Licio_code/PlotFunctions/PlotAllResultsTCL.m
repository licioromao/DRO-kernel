function [] = PlotAllResultsTCL(TCL_pathFiles,indexPaths,MC)

if isempty(indexPaths)
    L = length(TCL_pathFiles);
else
    L = indexPaths;
end

for i=1:L
    tempPlots = MonteCarloSimulation(TCL_pathFiles,'TCL',MC); 
    PlotResults(tempPlots);
end

[Data,EpValues] = GetKernelEpsilonsAndValueFunc(TCL_pathFiles,indexPaths);

L = length(Data);

if ~isempty(Data)
    figure
    hold on
    for i =1:L
        plot(Data{i}.grid_x,Data{i}.ValueFunc,'LineWidth',2);
    end
    box on
    legend(EpValues)
end

end

function [out,LegendNames] = GetKernelEpsilonsAndValueFunc(TCL_pathFiles,indexPaths)

% This function organises all the results and output the Kernel value
% function for different values of ep. If no Kernel ambiguity set is found
% in the given files the output is the empty string.
 
if isempty(indexPaths)
    L = length(TCL_pathFiles);
else
    L = indexPaths;
end

varnames = {'ValueFuncKernel','Grid'};
out = [];
LegendNames = [];

for i =1:L
    S{i} = load(TCL_pathFiles(i).FullPath,varnames{:});

    if ~isfield(S{i},'ValueFuncKernel') % Testing if there is ValueFuncKernel in the file
        warning('No ValueFuncKernel variable found in %s',TCL_pathFiles(i).FullPath);
    else % in the positive case, save the parameters
        out{i}.ep = S{i}.ValueFuncKernel.ParamAmbiguity.ep;
        LegendNames{i} = sprintf('%.5f',S{i}.ValueFuncKernel.ParamAmbiguity.ep);
        out{i}.gamma = S{i}.ValueFuncKernel.ParamAmbiguity.gamma;
        out{i}.ValueFunc =  S{i}.ValueFuncKernel.ValueFunctionQP(:,1);
        out{i}.grid_x = S{i}.Grid.getValues.Partition.grid_x;
    end 
end


end

