function out = getMomentPolicyTCL(horizon,numberPointsPartition ...
     ,numberDataEstimationTransition,ep,rhoMu,rhoSigma)
 
pathProject = getProjectPath;

temp = testTCLFunc(horizon,numberPointsPartition,numberDataEstimationTransition,ep,rhoMu,rhoSigma,pathProject);

load(temp,'TCL_ResultsPath')
load(TCL_ResultsPath(end).FullPath,'ValueFuncMoment','Grid')

out.OptPolicy = ValueFuncMoment.OptInput(:,1);
out.ValueFunc = ValueFuncMoment.ValueFunction(:,1);
out.horizon = horizon;
out.StatePartExp = Grid;
out.grid_x = Grid.getValues.Partition.grid_x;
out.param = ValueFuncMoment.param;
out.SafeSet = ValueFuncMoment.param.SafeSet;


end