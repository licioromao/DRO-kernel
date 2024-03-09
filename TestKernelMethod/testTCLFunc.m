function out = testTCLFunc(horizon,TCLpartition,mTCL,ep,rhoMu,rhoSigma,pathProject)

N = int16(horizon); 

% Empty to enable printing information of the script
StructNoAmbiguity.Name = [];
StructKernelAmbiguity.Name = [];
StructMomentAmbiguity.Name = [];
StructKLdivAmbiguity.Name = [];


if ~isempty(mTCL)
    N_TCL = size(mTCL,2);
    total_iterations = N_TCL;
else
    N_TCL = 0;
    total_iterations = 1;
end


if ~isempty(ep)
    N_ep = size(ep,2);
    total_iterations = total_iterations*N_ep;
else
    N_ep = 0;
end


if ~isempty(rhoMu)
    N_Mu = size(rhoMu,2);
    total_iterations = total_iterations*N_Mu;
else
    N_Mu = 0;
end


if ~isempty(rhoSigma)
    N_Sigma = size(rhoSigma,2);
    total_iterations = total_iterations*N_Sigma;
else
    N_Sigma = 0;
end

TCL_ResultsPath = cell(total_iterations,1);

% In addition to not passing as parameter in the function below, you should
% also comment here to ommit any ambiguity set
StructNoAmbiguity.Name = 'NoAmbiguity';
StructKernelAmbiguity.Name = 'KernelAmbiguity';
StructMomentAmbiguity.Name = 'MomentAmbiguity';
%StructKLdivAmbiguity.Name = 'KLdivAmbiguity';


timeIteration = 2; % A guess on the first iteration
OuterLoopInfo.TypeVectorField = 'TCL';
OuterLoopInfo.pathProject = pathProject;


StructKernelAmbiguity.gamma = 10;

OuterLoopInfo.TotalIteration = total_iterations;
OuterLoopInfo.TimeIteration = timeIteration;
TCL_ResultsPath = [];


for i1=1:N_TCL
    param = ComputeTransition(OuterLoopInfo.TypeVectorField,TCLpartition,mTCL(i1));
    index = RemainingIterations(1,[i1,N_TCL],total_iterations/N_TCL,[]);

    OuterLoopInfo.CurrentIteration = index;
    temp = TCL(N,{StructNoAmbiguity},OuterLoopInfo,param); % Solve the DP iteration without ambiguity set
    TCL_ResultsPath = [TCL_ResultsPath;temp];

    if ~isempty(ep)
        StructKernelAmbiguity.m = mTCL(i1);
        for i2=1:N_ep
            % Parameters of the kernel ambiguity set
            timeIterarionKey = tic;
            StructKernelAmbiguity.ep = ep(i2);
            index = RemainingIterations(2,[[i1;i2],[N_TCL;N_ep]],N_Mu*N_Sigma,[]);
            OuterLoopInfo.CurrentIteration = index;

            temp = TCL(N,{StructKernelAmbiguity},OuterLoopInfo,param); % Solving DP with Kernel Ambiguity set
            TCL_ResultsPath = [TCL_ResultsPath;temp];

            timeIteration = toc(timeIterarionKey); % Estimating how long time it took the last iteration
        end
    end

    if ~isempty(rhoMu) && ~isempty(rhoSigma)

        for i3=1:N_Mu
            for i4=1:N_Sigma


                % Parameters of the moment ambiguity set
                StructMomentAmbiguity.rhoMu = rhoMu(i3);
                StructMomentAmbiguity.rhoSigma = rhoSigma(i4);

                % Parameters of the kernel ambiguity set
                timeIterarionKey = tic;

                index = RemainingIterations(2,[[i3;i4],[N_Mu;N_Sigma]],1,[]);
                OuterLoopInfo.CurrentIteration = index;

                temp = TCL(N,{StructMomentAmbiguity},OuterLoopInfo,param); % Solving DP problem with Moment ambiguity set
                TCL_ResultsPath = [TCL_ResultsPath;temp];

                timeIteration = toc(timeIterarionKey); % Estimating how long time it took the last iteration

            end
        end
    end
end


out = strcat(pathProject, ... 
    sprintf('/Results/results_%s/TCL/paths_%s.mat',char(java.net.InetAddress.getLocalHost.getHostName), ... 
    TCL_ResultsPath(end).FileName));

save(out)

end
