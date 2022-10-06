function out = MonteCarloSimulation(PathToFile,TypeOfVectorField,MC)


switch TypeOfVectorField % checking the type of ambiguity
    case 'TCL'
        index = strfind(PathToFile,'TCL');
        if isempty(index)
            error('The path you provided does not TCL data'); % error if TCL is not present in the path
        end
        
        ValueFunc = load(PathToFile,'-regexp','ValueFunc*');
        ValueFunc = struct2cell(ValueFunc);
        L = length(ValueFunc);
        
        load(PathToFile,'StructAmbiguityTypes');
        
        for i=1:L
            if isfield(StructAmbiguityTypes{i},'ep')
                AmbiguityParam.ep = StructAmbiguityTypes{i}.ep;
            end
            
            if isfield(StructAmbiguityTypes{i},'rhoMu')
                AmbiguityParam.rhoMu = StructAmbiguityTypes{i}.rhoMu;
            end
            
            if isfield(StructAmbiguityTypes{i},'rhoSigma')
                AmbiguityParam.rhoSigma = StructAmbiguityTypes{i}.rhoSigma;
            end
        end
        
        load(PathToFile,'param')
        
        m = param.MC;
        PartitionSize = param.NumberOfPartitions;
         
        if L == 0
            error('There is no value function on this data file'); % error if there is not value function 
        end
        
        Partition = load(PathToFile,'Grid');
        load(PathToFile,'TimeHorizon');
        load(PathToFile,'param');
        
        Partition = Partition.Grid;
        grid_x = Partition.getValues.Partition.grid_x;
        
        NumberOfPoints = size(grid_x,1);     
        
        tic
        for i = 1:L
            
            sum_Safety = zeros(NumberOfPoints,1);
            Inputs = ValueFunc{i}.OptInput;
            
            if strcmp(ValueFunc{i}.AmbiguityType,'KernelAmbiguity')
                sum_Safety_2 = zeros(NumberOfPoints,1);
                InputsConservative = ValueFunc{i}.OptInputConservative;
            end
            
            
            for k=1:NumberOfPoints
                x0 = grid_x(k);
                for j = 1:MC
                    sum_Safety(k) = sum_Safety(k) + SimulateTrajectory(x0,Partition,Inputs,TimeHorizon,'TCL',param);
                    sum_Safety_2(k) = sum_Safety_2(k) + SimulateTrajectory(x0,Partition,InputsConservative,TimeHorizon,'TCL',param);
                end
            end
            
            if strcmp(ValueFunc{i}.AmbiguityType,'KernelAmbiguity')
                out{i}.ValueFunctionConservative = ValueFunc{i}.ValueFunctionConservative(:,1);
                out{i}.OptInputConservative = ValueFunc{i}.OptInputConservative;
                out{i}.EmpiricalValueFuncConservative = sum_Safety_2/MC;
            end
            
            out{i}.ValueFunction = ValueFunc{i}.ValueFunction(:,1);
            out{i}.OptInput = ValueFunc{i}.OptInput;
            out{i}.AmbiguityType = ValueFunc{i}.AmbiguityType;
            out{i}.TypeOfVectorField = ValueFunc{i}.TypeOfVectorField;
            out{i}.TimeHorizon = ValueFunc{i}.N;
            out{i}.EmpiricalValueFunc = sum_Safety/MC;
            
        end
            toc
        
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

