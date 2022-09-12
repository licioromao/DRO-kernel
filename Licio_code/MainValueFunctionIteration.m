% Function that computes the value function

function ValueFunc = MainValueFunctionIteration(PartitionObj,InputPartition,TypeOfVectorField,TypeOfAmbiguity,BooleanSignal,param)

if ~BooleanSignal
    switch TypeOfAmbiguity
        case 'NoAmbiguity'
            fprintf('\nComputing value function (without ambiguity)...\n')
        case 'MomentAmbiguity'
            fprintf('\nComputing value function (moment ambiguity)...\n')
        case 'KernelAmbiguity'
            fprintf('\nComputing value function (kernel ambiguity)...\n')
        case 'KLdivAmbiguity'
            fprintf('\nComputing value function (KL ambiguity)...\n')
        case 'WassersteinAmbiguity'
            fprintf('\nComputing value function (Wasserstein ambiguity)...\n')
        otherwise
            error('Type of Ambiguity not implemented')
    end
    
    tic;
    ValueFunc = ComputeValueFunction(param,TypeOfVectorField,TypeOfAmbiguity);
    
    ValueFunc = ValueFunc.getIndexReachAvoid(PartitionObj.getValues);
    ValueFunc = ValueFunc.BackwardIteration(PartitionObj,InputPartition);
    ValueFunc.time = toc;
    
    fprintf('Done\n');
else
    ValueFunc = [];
    switch TypeOfAmbiguity
        case 'NoAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncNoAmbiguity');
        case 'MomentAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncMoment');
        case 'KernelAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncKernel');
        case 'KLdivAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncKLdivAmbiguity');
        case 'WassersteinAmbiguity'
            warning('The variable %s already exists. Skipping to the next string...\n','ValueFuncWasserstein');
        otherwise
            error('Type of Ambiguity not implemented')
    end
end








end