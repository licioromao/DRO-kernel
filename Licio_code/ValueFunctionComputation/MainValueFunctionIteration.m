% Function that computes the value function

function ValueFunc = MainValueFunctionIteration(PartitionObj,InputPartition,TypeOfVectorField,StructAmbiguityTypes,BooleanSignal,param)

if ~BooleanSignal
    switch StructAmbiguityTypes.Name
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
    
    ValueFunc = ComputeValueFunction(param,TypeOfVectorField,StructAmbiguityTypes);
    
    switch TypeOfVectorField
        case 'Fishery'
            ValueFunc = ValueFunc.getIndexReachAvoid(PartitionObj.getValues);
        case 'TCL'
            ValueFunc = ValueFunc.getIndexSafety(PartitionObj.getValues);
        case 'ChainInt'
            ValueFunc = ValueFunc.getIndexSafety(PartitionObj.getValues);
        case 'CarPole'
            ValueFunc = ValueFunc.getIndexReachAvoid(PartitionObj.getValues);
        case 'CarPoleNL'
            ValueFunc = ValueFunc.getIndexReachAvoid(PartitionObj.getValues);
        otherwise
            NotImplemented();
    end
    
    ValueFunc = ValueFunc.BackwardIteration(PartitionObj,InputPartition);
    
    fprintf('Done\n');
else
    ValueFunc = [];
    switch StructAmbiguityTypes.Name
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