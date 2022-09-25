function [] = PlotAllResultsTCL(TCL_pathFiles,indexPaths,MC)

if isempty(indexPaths)
    L = length(TCL_pathFiles);
else
    L = indexPaths;
end

for i=1:L
    tempPlots = MonteCarloSimulation(TCL_pathFiles{i}.FullPath,'TCL',MC); 
    PlotResults(tempPlots);
end


end

