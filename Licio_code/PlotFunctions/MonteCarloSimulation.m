function out = MonteCarloSimulation(PathToFile,TypeOfVectorField,MC)

% This function collects all the data and plot the corresponding value
% function

switch TypeOfVectorField % checking the type of ambiguity
    case 'TCL'
        index = strfind(PathToFile(1).FullPath,'TCL');
        if isempty(index)
            error('The path you provided does not TCL data'); % error if TCL is not present in the path
        end
        
        ValueFunc = [];
        StructAmbiguityTypes = [];
        for i = 1:length(PathToFile)
            ValueFunc{i} = load(PathToFile(i).FullPath,'-regexp','ValueFunc*');
            StructAmbiguityTypes{i} = load(PathToFile(i).FullPath,'StructAmbiguityTypes');
        end
        
        L = length(ValueFunc);

        if L == 0
            error('There is no value function on this data file'); % error if there is not value function 
        end

        indexKernel = 0;
        indexMoment = 0;
        
        for i=1:L
            L2 = length(StructAmbiguityTypes{i});
            for j=1:L2
                if isfield(StructAmbiguityTypes{i}.StructAmbiguityTypes{j},'ep')
                    indexKernel = indexKernel + 1;
                    AmbiguityParam.ep(indexKernel) = StructAmbiguityTypes{i}.StructAmbiguityTypes{j}.ep;
                end

                if isfield(StructAmbiguityTypes{i}.StructAmbiguityTypes{j},'rhoMu')
                    indexMoment = indexMoment + 1;
                    AmbiguityParam.rhoMu(indexMoment) = StructAmbiguityTypes{i}.StructAmbiguityTypes{j}.rhoMu;
                    AmbiguityParam.rhoSigma(indexMoment) = StructAmbiguityTypes{i}.StructAmbiguityTypes{j}.rhoSigma;
                end
            end
        end
        
        load(PathToFile(1).FullPath,'param')
        
        m = param.MC;
        PartitionSize = param.NumberOfPartitions;
        
        Partition = load(PathToFile(1).FullPath,'Grid');
        load(PathToFile(1).FullPath,'TimeHorizon');
        load(PathToFile(1).FullPath,'param');
        
        Partition = Partition.Grid;
        grid_x = Partition.getValues.Partition.grid_x;
        
        NumberOfPoints = size(grid_x,1);     
        
        t_start = tic;
        for i = 1:L
            L2 = length(fieldnames(ValueFunc{i}));
            NamesField = fieldnames(ValueFunc{i});

            for j = 1:L2
                index = RemainingIterations(2,[[i;j],[L;L2]],1,[]);
                ParamAmbiguity = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.ParamAmbiguity'));

                if strcmp(ParamAmbiguity.Name,'KernelAmbiguity')
                    Inputs{index} = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.OptInput'));
                    InputsConservative{index} = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.OptInputConservative'));
                    InputsQP{index} = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.OptInputQP'));
                    
                    Safety_2 = zeros(NumberOfPoints,1);
                    Safety_3 = zeros(NumberOfPoints,1);
                else
                    Inputs{index} = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.OptInput'));
                end 
                InfoValueFunc{index} = ParamAmbiguity;
                Safety = zeros(NumberOfPoints,1);

                for k1=1:NumberOfPoints
                    x0 = grid_x(k1);

                    if strcmp(InfoValueFunc{index}.Name,'KernelAmbiguity')
                        [Safety(k1),Safety_2(k1),Safety_3(k1)] = ...
                            ComputeEmpiricalProb(x0,Partition,Inputs{index},InputsConservative{index},InputsQP{index}...
                            ,TimeHorizon,InfoValueFunc{index}.Name, ...
                            'TCL',MC,param);
                    else
                        [Safety(k1),~,~] = ...
                            ComputeEmpiricalProb(x0,Partition,Inputs{index},[],[],TimeHorizon,InfoValueFunc{index}.Name, ...
                            'TCL',MC,param);
                    end
                end

                out{index}.ValueFunction = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.ValueFunction(:,1)'));
                out{index}.OptInput = Inputs{index};
                out{index}.EmpiricalValueFunc = Safety;

                if strcmp(InfoValueFunc{index}.Name,'KernelAmbiguity')
                    out{index}.ValueFunctionConservative = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.ValueFunctionConservative(:,1)'));
                    out{index}.OptInputConservative = InputsConservative{index};
                    out{index}.EmpiricalValueFuncConservative = Safety_2;

                    out{index}.ValueFunctionMatrix = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.ValueFunction(:,1)'));
                    out{index}.OptInputMatrix = Inputs{index};
                    out{index}.EmpiricalValueFuncMatrix = Safety;

                    out{index}.ValueFunction = eval(strcat(sprintf('ValueFunc{%d}.',i),NamesField{j},'.ValueFunctionQP(:,1)'));
                    out{index}.OptInput = InputsQP{index};
                    out{index}.EmpiricalValueFunc = Safety_3;
                end

                out{index}.AmbiguityType = InfoValueFunc{index}.Name;
                out{index}.TypeOfVectorField = 'TCL';
                out{index}.TimeHorizon = TimeHorizon;
            end
             
        end
        toc(t_start)
        
    otherwise
        NotImplemented(); % error if vector field is not implemented
end


out{L+1}.Partition = Partition;
out{L+1}.AmbiguityParam = AmbiguityParam;
out{L+1}.m = m;
out{L+1}.PartitionSize = PartitionSize;

end

function out = SimulateTrajectory(InitialState,Partition,Inputs,TimeHorizon,TypeOfVectorField,param)

switch TypeOfVectorField
    case 'TCL'
        SafeSet = param.SafeSet;
        out = 1;
        
        CurrentState = InitialState;
        
        tempPartition = Partition.getElementPartition(CurrentState);
        u = Inputs(tempPartition.index);
        
        VectorFieldObj = VectorFieldTCL(CurrentState,u,param);
        
        if CurrentState < SafeSet(1) || CurrentState > SafeSet(2)
            out = 0;
            return;
        else
            for k = 1:TimeHorizon
                VectorFieldObj.Noise = random(param.w);
                CurrentState = VectorFieldObj.IterateDynamics;
                
                if CurrentState < SafeSet(1) || CurrentState > SafeSet(2)
                    out = 0;
                    return;
                end
                
                tempPartition = Partition.getElementPartition(CurrentState);
                u = Inputs(tempPartition.index);
                
                VectorFieldObj.SetValues(CurrentState,u,param);
            end
        end
        
    otherwise
        NotImplemented();
end

end

function [sum_Safety,sum_Safety_2,sum_Safety_3] = ComputeEmpiricalProb(x,Partition,OptInput,OptInputConservative...
     ,OptInputQP,TimeHorizon,AmbiguityType,VectorField,MC,param)

sum_Safety = 0;

if strcmp(AmbiguityType,'KernelAmbiguity')
    sum_Safety_2 = 0;
    sum_Safety_3 = 0;
end

if strcmp(AmbiguityType,'KernelAmbiguity')
    for k = 1:MC
        sum_Safety = sum_Safety + SimulateTrajectory(x,Partition,OptInput,TimeHorizon,VectorField,param);
        sum_Safety_2 = sum_Safety_2 + SimulateTrajectory(x,Partition,OptInputConservative,TimeHorizon,VectorField,param);
        sum_Safety_3 = sum_Safety_3 + SimulateTrajectory(x,Partition,OptInputQP,TimeHorizon,VectorField,param);
    end
else
    for k = 1:MC
        sum_Safety = sum_Safety + SimulateTrajectory(x,Partition,OptInput,TimeHorizon,VectorField,param);
    end
end

sum_Safety = sum_Safety/MC;

if strcmp(AmbiguityType,'KernelAmbiguity')
    sum_Safety_2 = sum_Safety_2/MC;
    sum_Safety_3 = sum_Safety_3/MC;
else
    sum_Safety_2 = [];
    sum_Safety_3 = [];
end

end
