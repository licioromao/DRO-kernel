% Empty to enable printing information of the script
StructNoAmbiguity.Name = [];
StructKernelAmbiguity.Name = [];
StructMomentAmbiguity.Name = [];
StructKLdivAmbiguity.Name = [];

% Time horizon
N = int16(8);

% Number of points between 18 and 24 degree 
TCLpartition = 2000;

mTCL = 5000;
ep = [0.001,0.01,0.1,1];
rhoMu = [0.1,0.5,0.7];
rhoSigma = [0.5,1,2,5];

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


timeIteration = 40; % A guess on the first iteration
OuterLoopInfo.TypeVectorField = 'TCL';

for i1=1:N_TCL
    for i2=1:N_ep
        for i3=1:N_Mu
            for i4=1:N_Sigma
                
                timeIterarionKey = tic;
                index = RemainingIterations(4,[[i1;i2;i3;i4],[N_TCL;N_ep;N_Mu;N_Sigma]],1,[]);
                
                OuterLoopInfo.TotalIteration = total_iterations;
                OuterLoopInfo.CurrentIteration = index;
                OuterLoopInfo.TimeIteration = timeIteration;
                
                % Parameters of the kernel ambiguity set
                StructKernelAmbiguity.ep = ep(i2);
                StructKernelAmbiguity.gamma = 1;
                StructKernelAmbiguity.m = mTCL(i1);
                
                % Parameters of the moment ambiguity set            
                StructMomentAmbiguity.rhoMu = rhoMu(i3);
                StructMomentAmbiguity.rhoSigma = rhoSigma(i4);
                
                % Parameters of the KL ambiguity set
                StructKLdivAmbiguity.ep = ep(i2);               
                
                TCL_ResultsPath{index} = TCL(N,TCLpartition,mTCL(i1),{StructNoAmbiguity,StructKernelAmbiguity,StructMomentAmbiguity,StructKLdivAmbiguity},OuterLoopInfo);
                
                timeIteration = toc(timeIterarionKey); % Estimating how long time it took the last iteration
                
                %TCL_ResultsPath{index} = TCL(N,TCLpartition,mTCL(i1),{StructNoAmbiguity,StructMomentAmbiguity,StructKernelAmbiguity});
                
            end
        end
    end
end







