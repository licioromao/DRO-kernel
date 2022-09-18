factor = 1000;

% Empty to enable printing information of the script
StructNoAmbiguity.Name = [];
StructKernelAmbiguity.Name = [];
StructMomentAmbiguity.Name = [];
StructKLdivAmbiguity.Name = [];

% Time horizon
N = int16(2);

% Partition of the state space
FisheryPartition = [2,2,2];

mFishery = [200];
ep = [0.1];
rhoMu = [10];
rhoSigma = [5];

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

% % Creating a progress bar of the value function computation
% hh = waitbar(0,'Outer loop iterations','Name','Outer loop');
% total_iterations = N_Fishery*N_ep*N_Mu*N_Sigma;
% fprintf('(Outer loop) Total of iterations: %d \n',total_iterations);

timeIteration = 40; % A guess on the first iteration
OuterLoopInfo.TypeVectorField = 'Fishery';

for i1=1:N_Fishery
    for i2=1:N_ep
        for i3=1:N_Mu
            for i4=1:N_Sigma
                
                timeIterarionKey = tic;
                index = RemainingIterations(4,[[i1;i2;i3;i4],[N_TCL;N_ep;N_Mu;N_Sigma]],1,[]);
                
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
               
                
                %Fishery_ResultsPath{index} = Fishery_management(N,FisheryPartition,mFishery(i1),{StructKLdivAmbiguity,StructNoAmbiguity,StructMomentAmbiguity,StructKernelAmbiguity});
                %Fishery_ResultsPath{index} = Fishery_management(N,FisheryPartition,mFishery(i1),{StructNoAmbiguity,StructKernelAmbiguity});
                Fishery_ResultsPath{index} = Fishery_management(N,FisheryPartition,mFishery(i1),{StructNoAmbiguity,StructMomentAmbiguity},OuterLoopInfo);

%                 IterationsToGo = (total_iterations - index);
%                 perc_iterates = index/total_iterations;
%                 waitbar(perc_iterates,hh,sprintf('%.5f completed. %d iterations to go',perc_iterates,IterationsToGo));
%                 
                
            end
        end
    end
end

