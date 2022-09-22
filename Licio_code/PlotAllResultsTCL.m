function [] = PlotAllResultsTCL(TCL_pathFiles,MC)

L = length(TCL_pathFiles);

for i=1:L
    tempPlots = MonteCarloSimulation(TCL_pathFiles{i}.FullPath,'TCL',MC); 
    PlotResults(tempPlots);
end


end

