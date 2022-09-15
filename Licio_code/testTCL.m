factor = 1000;

N = int16(18);

TCLpartition = 2;

mTCL = [5,15];
ep = 0.01;
rhoMu = [0.1];
rhoSigma = [1,2];

N_TCL = size(mTCL,2);
N_ep = size(ep,2);
N_Mu = size(rhoMu,2);
N_Sigma = size(rhoSigma,2);

TCL_ResultsPath = cell(N_TCL*N_ep*N_Mu*N_Sigma,1);

StructNoAmbiguity.Name = 'NoAmbiguity';


% Creating a progress bar of the value function computation
hh = waitbar(0,'Outer loop iterations','Name','Outer loop');
total_iterations = N_TCL*N_ep*N_Mu*N_Sigma;
fprintf('(Outer loop) Total of iterations: %d \n',total_iterations);

for i1=1:N_TCL
    for i2=1:N_ep
        for i3=1:N_Mu
            for i4=1:N_Sigma
                
                % Structure with the Kernel
                StructKernelAmbiguity.Name = 'KernelAmbiguity';
                StructKernelAmbiguity.ep = ep(i2);
                StructKernelAmbiguity.gamma = 1;
                StructKernelAmbiguity.m = mTCL(i1);
                
                % Structure with Moment Ambiguity sets                
                StructMomentAmbiguity.Name = 'MomentAmbiguity';
                StructMomentAmbiguity.rhoMu = rhoMu(i3);
                StructMomentAmbiguity.rhoSigma = rhoSigma(i4);
                
                
                StructKLdivAmbiguity.Name = 'KLdivAmbiguity';
                StructKLdivAmbiguity.ep = ep(i2);
                
                index = RemainingIterations(4,[[i1;i2;i3;i4],[N_TCL;N_ep;N_Mu;N_Sigma]],1,hh);
                
                TCL_ResultsPath{index} = TCL(N,TCLpartition,mTCL(i1),{StructKLdivAmbiguity,StructNoAmbiguity,StructMomentAmbiguity,StructKernelAmbiguity});
                
                IterationsToGo = (total_iterations - index);
                perc_iterates = index/total_iterations;
                waitbar(perc_iterates,hh,sprintf('%.5f completed. %d iterations to go',perc_iterates,IterationsToGo));
                
            end
        end
    end
end

delete(hh)







