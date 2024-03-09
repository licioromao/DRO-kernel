AddFunctionPaths();

% Empty to enable printing information of the script
StructNoAmbiguity.Name = [];
StructKernelAmbiguity.Name = [];
StructMomentAmbiguity.Name = [];
StructKLdivAmbiguity.Name = [];

% Time horizon
N = int16(10);

% Partition of the state space
FisheryPartition = [2,2,2];

mFishery = [110];
ep = [0.1];
rhoMu = [300];
rhoSigma = [50];

N_Fishery = size(mFishery,2);
N_ep = size(ep,2);
N_Mu = size(rhoMu,2);
N_Sigma = size(rhoSigma,2);

total_iterations = N_Fishery*N_ep*N_Mu*N_Sigma;

Fishery_ResultsPath = cell(total_iterations,1);

StructNoAmbiguity.Name = 'NoAmbiguity';
StructKernelAmbiguity.Name = 'KernelAmbiguity';
StructMomentAmbiguity.Name = 'MomentAmbiguity';
StructKLdivAmbiguity.Name = 'KLdivAmbiguity';

timeIteration = 40; % A guess on the first iteration
OuterLoopInfo.TypeVectorField = 'Fishery';

for i1=1:N_Fishery
    for i2=1:N_ep
        for i3=1:N_Mu
            for i4=1:N_Sigma
                
                timeIterarionKey = tic;
                index = RemainingIterations(4,[[i1;i2;i3;i4],[N_Fishery;N_ep;N_Mu;N_Sigma]],1,[]);
                
                OuterLoopInfo.TotalIteration = total_iterations;
                OuterLoopInfo.CurrentIteration = index;
                OuterLoopInfo.TimeIteration = timeIteration;
                
                % Structure with the Kernel
                StructKernelAmbiguity.Name = 'KernelAmbiguity';
                StructKernelAmbiguity.ep = ep(i2);
                StructKernelAmbiguity.gamma = 1;
                StructKernelAmbiguity.m = mFishery(i1);
                
                % Structure with Moment Ambiguity sets                
                StructMomentAmbiguity.Name = 'MomentAmbiguity';
                StructMomentAmbiguity.rhoMu = rhoMu(i3);
                StructMomentAmbiguity.rhoSigma = rhoSigma(i4);
                
                
                StructKLdivAmbiguity.Name = 'KLdivAmbiguity';
                StructKLdivAmbiguity.ep = ep(i2);
               
                
                Fishery_ResultsPath{index} = Fishery_management(N,FisheryPartition,mFishery(i1),{StructKLdivAmbiguity,StructNoAmbiguity,StructMomentAmbiguity,StructKernelAmbiguity},OuterLoopInfo);
                %Fishery_ResultsPath{index} = Fishery_management(N,FisheryPartition,mFishery(i1),{StructNoAmbiguity,StructKernelAmbiguity});
%                Fishery_ResultsPath{index} = Fishery_management(N,FisheryPartition,mFishery(i1),{StructNoAmbiguity},OuterLoopInfo);
                 
                
            end
        end
    end
end

PathsString = sprintf('./Results/results_%s/Fishery/paths_%s.mat',char(java.net.InetAddress.getLocalHost.getHostName),TCL_ResultsPath{end}.FileName);

save(PathsString);

RemoveFunctionPaths();
