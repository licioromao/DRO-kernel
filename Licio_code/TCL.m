function out = TCL(TimeHorizon,StructAmbiguityTypes,OuterLoopInfo,InputParam)

param = InputParam;

param.OuterLoopInfo = OuterLoopInfo;
param.OuterLoopInfo.StringAmbiguity = StructAmbiguityTypes; % THIS IS BAD PRACTICE SINCE I AM ALSO SHARING AMBIGUITY PARAMETERS
param.N = int8(TimeHorizon);

Grid = StatePartition(param.NumberOfPartitions,param.SafeSet,'TCL'); % Generate the partition of the state space
InputPartition = generateInputPartition([],'TCL'); % Generate a vector with all possible combinations of inputs


L = length(StructAmbiguityTypes);

for i = 1:L
    switch StructAmbiguityTypes{i}.Name
        case 'NoAmbiguity'
            TotalTime = tic;
            ValueFuncNoAmbiguity = MainValueFunctionIteration(Grid,InputPartition,'TCL',StructAmbiguityTypes{i},exist('ValueFuncNoAmbiguity','var'),param);
            ValueFuncNoAmbiguity.time = toc(TotalTime);
            paramSave = [];
            
        case 'MomentAmbiguity'
            TotalTime1 = tic;
            ValueFuncMoment = MainValueFunctionIteration(Grid,InputPartition,'TCL',StructAmbiguityTypes{i},exist('ValueFuncMoment','var'),param);
            ValueFuncMoment.time = toc(TotalTime1);
            
            paramSave.rhoMu = StructAmbiguityTypes{i}.rhoMu;
            paramSave.rhoSigma = StructAmbiguityTypes{i}.rhoSigma;
             
        case 'WassersteinAmbiguity'
            TotalTime2 = tic;
            ValueFuncWasserstein = MainValueFunctionIteration(Grid,InputPartition,'TCL',StructAmbiguityTypes{i},exist('ValueFuncWasserstein','var'),param);
            ValueFuncWasserstein.time = toc(TotalTime2);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
                   
            
        case 'KLdivAmbiguity'
            TotalTime3 = tic;
            ValueFuncKL = MainValueFunctionIteration(Grid,InputPartition,'TCL',StructAmbiguityTypes{i},exist('ValueFuncKL','var'),param);
            ValueFuncKL.time = toc(TotalTime3);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
            
        case 'KernelAmbiguity'  
            TotalTime4 = tic;
            ValueFuncKernel = MainValueFunctionIteration(Grid,InputPartition,'TCL',StructAmbiguityTypes{i},exist('ValueFuncKernel','var'),param);
            ValueFuncKernel.time = toc(TotalTime4);
            
            paramSave.ep = StructAmbiguityTypes{i}.ep;
           
        otherwise
            warning('%s has not been implemented. Jumping to the next string',StructAmbiguityTypes{i});
    end   
end

FileName = getDateSaveFile(TimeHorizon,param.NumberOfPartitions,param.MC,'TCL',paramSave); % Getting the name of file based on the current date and time
save(FileName.FullPath); % saving the results in the path specified by FILE

out = FileName;


end

