clear all
close all

AddFunctionPaths;

% Empty to enable printing information of the script
StructNoAmbiguity.Name = [];
StructKernelAmbiguity.Name = [];
StructMomentAmbiguity.Name = [];
StructKLdivAmbiguity.Name = [];

% Time horizon
N = int16(3);

% Number of points between 18 and 24 degree 
TCLpartition = 30;

mTCL = 100;
% 0.01, 0.05, 0.5
ep = [0,1e-4,0.01];
rhoMu = 0.7;
rhoSigma = [3];

N_TCL = size(mTCL,2);
N_ep = size(ep,2);
N_Mu = size(rhoMu,2);
N_Sigma = size(rhoSigma,2);

total_iterations = N_TCL*N_ep*N_Mu*N_Sigma;

TCL_ResultsPath = cell(total_iterations,1);

% Inaddition to not passing as parameter in the function below, you should
% also comment here to ommit any ambiguity set
StructNoAmbiguity.Name = 'NoAmbiguity';
StructKernelAmbiguity.Name = 'KernelAmbiguity';
StructMomentAmbiguity.Name = 'MomentAmbiguity';
StructKLdivAmbiguity.Name = 'KLdivAmbiguity';


timeIteration = 2; % A guess on the first iteration
OuterLoopInfo.TypeVectorField = 'TCL';


StructKernelAmbiguity.gamma = 10;

OuterLoopInfo.TotalIteration = total_iterations;
OuterLoopInfo.TimeIteration = timeIteration;
TCL_ResultsPath = [];


for i1=1:N_TCL
    param = ComputeTransition(OuterLoopInfo.TypeVectorField,TCLpartition,mTCL(i1));
    index = RemainingIterations(1,[i1,N_TCL],N_Mu*N_Sigma*N_ep,[]);

    OuterLoopInfo.CurrentIteration = index;
    temp = TCL(N,{StructNoAmbiguity},OuterLoopInfo,param);
    StructKernelAmbiguity.m = mTCL(i1);
    for i2=1:N_ep
        % Parameters of the kernel ambiguity set
        timeIterarionKey = tic;
        StructKernelAmbiguity.ep = ep(i2);
        index = RemainingIterations(2,[[i1;i2],[N_TCL;N_ep]],N_Mu*N_Sigma,[]);
        OuterLoopInfo.CurrentIteration = index;

        temp = TCL(N,{StructKernelAmbiguity},OuterLoopInfo,param);
        TCL_ResultsPath = [TCL_ResultsPath;temp];

        timeIteration = toc(timeIterarionKey); % Estimating how long time it took the last iteration
    end

    for i3=1:N_Mu
        for i4=1:N_Sigma


            % Parameters of the moment ambiguity set
            StructMomentAmbiguity.rhoMu = rhoMu(i3);
            StructMomentAmbiguity.rhoSigma = rhoSigma(i4);

            % Parameters of the kernel ambiguity set
            timeIterarionKey = tic;
            
            index = RemainingIterations(2,[[i3;i4],[N_Mu;N_Sigma]],1,[]);
            OuterLoopInfo.CurrentIteration = index;

            temp = TCL(N,{StructMomentAmbiguity},OuterLoopInfo,param);
            TCL_ResultsPath = [TCL_ResultsPath;temp];

            timeIteration = toc(timeIterarionKey); % Estimating how long time it took the last iteration

        end
    end
end


PathsString = sprintf('./Results/results_%s/TCL/paths_%s.mat',char(java.net.InetAddress.getLocalHost.getHostName),TCL_ResultsPath(end).FileName);

save(PathsString);

RemoveFunctionPaths;







