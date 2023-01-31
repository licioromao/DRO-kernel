
addpath ../AuxiliaryFunctions/
pathProject = getProjectPath;

AddFunctionPaths(pathProject);

generateNewData = 1;

horizon = 4;
numberOfPartition = 300;
numberOfDataTransition = 10000;
rhoMu = 1;
rhoSigma = 5;

if generateNewData == 1
    momentData = getMomentPolicyTCL(horizon,numberOfPartition,numberOfDataTransition ... 
        ,[],rhoMu,rhoSigma); % generating optimal policy for moment-based ambiguity
    save('momentData.mat',"momentData");
else
    load('momentData.mat')
end

lambdaValues = 2;
sigmaValues = 0.01;
eta = 1;

numberOfSimulationData = 500;

N1 = length(numberOfSimulationData);
N2 = length(lambdaValues);
N3 = length(sigmaValues);
N4 = length(eta);

for ell=1:N1
    
    figure
    hold on
    ValueFuncFinal = zeros(numberOfSimulationData(ell),N2*N3*N4);

    [DataCME,StatePartCME] = generateDataTCL(momentData.StatePartExp, ...
        momentData.OptPolicy,numberOfSimulationData(ell),momentData.param); % simulation data of the dynamics using the policy obtained in the previous step

    labelEstimatedVF = @(lambda,sigma,eta) valueFunctionKernelEstimation(momentData.horizon,...
        DataCME,StatePartCME,momentData.SafeSet,lambda,sigma,eta);

    for i=1:N2
        for j=1:N3
            tempValueFunc = zeros(numberOfSimulationData(ell),N4);
            for k=1:N4
                iterationsToGo = RemainingIterations(4,[[ell;i;j;k],[N1;N2;N3;N4]],1,[]);
                timecount = tic;
                ValueFuncFinal(:,iterationsToGo) = labelEstimatedVF(lambdaValues(i),sigmaValues(j),eta(k));
                PrintLoopProgress(N1*N2*N3*N4,iterationsToGo,toc(timecount));   
            end
        end
    end


    plot(DataCME(:,2),ValueFuncFinal,'LineWidth',2);
    plot(momentData.grid_x,momentData.ValueFunc,'LineWidth',2);
    axis([18,23,0,1])
end

RemoveFunctionPaths(pathProject);
