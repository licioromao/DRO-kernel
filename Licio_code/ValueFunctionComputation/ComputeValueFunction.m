classdef ComputeValueFunction
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        ValueFunction % This value stores the optimal value function
        ValueFunctionConservative % This value stores the optimal value function
        ValueFunctionQP % This value stores the optimal value function

        TypeOfVectorField % This value stores the type of Vector field
        
        OptInput % This value stores the optimal action
        OptInputConservative % This value stores the optimal action
        OptInputQP % This value stores the optimal action

        AmbiguityType % Type of ambiguity set
        ParamAmbiguity % Stores the parameters of the Ambiguity type
        
        param % structure with all the parameters of the problem
        
        IndexSafeAndReachSet % Indices of the safe and reach set
        IndexSafeSet
        
        N % horizon of the reach-avoid property
        
        time % this will store the time to go through a full value function computation
    end
    
    methods
        function obj = ComputeValueFunction(ArgParam,TypeOfVectorField,StructAmbiguityTypes)
            
            switch TypeOfVectorField
                case 'TCL'
                    NumberOfPoints = ArgParam.NumberOfPartitions; % number of points of the value function
                case 'ChainInt'
                    NumberOfPoints = ArgParam.NumberOfPartitions(1)*ArgParam.NumberOfPartitions(2); % number of points of the value function    
                case 'Fishery'
                    NumberOfPoints = ArgParam.NumberOfPartitions(1)*ArgParam.NumberOfPartitions(2)*ArgParam.NumberOfPartitions(3); % number of points of the value function
                case 'CarPole'
                    NumberOfPoints = ArgParam.NumberOfPartitions(1)*ArgParam.NumberOfPartitions(2)*ArgParam.NumberOfPartitions(3)*ArgParam.NumberOfPartitions(4); % number of points of the value function
                case 'CarPoleNL'
                    NumberOfPoints = ArgParam.NumberOfPartitions(1)*ArgParam.NumberOfPartitions(2)*ArgParam.NumberOfPartitions(3)*ArgParam.NumberOfPartitions(4); % number of points of the value function
                otherwise
                    NotImplemented();
            end
            
            
            obj.ValueFunction = zeros(NumberOfPoints,ArgParam.N+1); % initializing the variable to store the value function
            obj.ValueFunctionConservative = zeros(NumberOfPoints,ArgParam.N+1);  % This value stores the optimal value function for the matrix factorization case if ambiguity set is equal to kernel
            obj.ValueFunctionQP = zeros(NumberOfPoints,ArgParam.N+1);  % This value stores the optimal value function for the QP case if ambiguity set is equal to kernel

            obj.OptInput = zeros(NumberOfPoints,ArgParam.N+1); % initializing the variable to store the optimal input
            obj.OptInputConservative = zeros(NumberOfPoints,ArgParam.N+1); % initializing the variable to store the optimal input for kernel ambiguity, which has three value functions
            obj.OptInputQP = zeros(NumberOfPoints,ArgParam.N+1); % initializing the variable to store the optimal input for kernel ambiguity, which has three value functions
            
            obj.TypeOfVectorField = TypeOfVectorField;
            obj.AmbiguityType = StructAmbiguityTypes.Name;
            obj.ParamAmbiguity = StructAmbiguityTypes;
            
            obj.param = ArgParam;
            
            obj.IndexSafeAndReachSet = [];
            
            obj.N = ArgParam.N;
            
            obj.time = [];
        end
        
        function obj = getIndexReachAvoid(obj,ObjStatePartition)
            
            if ~(strcmp(obj.TypeOfVectorField,'Fishery') || strcmp(obj.TypeOfVectorField,'CarPole') || strcmp(obj.TypeOfVectorField,'CarPoleNL') )
                error('This function can only be used for reach avoid specifications. Please check if this is your case.')
            end
            
            % Retuns the indices of the safe and reach sets associated with the partition
            
            tempStatePartition = ObjStatePartition.Partition;
            
            switch obj.TypeOfVectorField
                case 'Fishery'
                    [obj.IndexSafeAndReachSet.safeIndex,obj.IndexSafeAndReachSet.reachIndex] = GetReachAvoid('3D',tempStatePartition,obj.param);
                case 'CarPole'
                    [obj.IndexSafeAndReachSet.safeIndex,obj.IndexSafeAndReachSet.reachIndex] = GetReachAvoid('4D',tempStatePartition,obj.param);
                case 'CarPoleNL'
                    [obj.IndexSafeAndReachSet.safeIndex,obj.IndexSafeAndReachSet.reachIndex] = GetReachAvoid('4D',tempStatePartition,obj.param);
                otherwise
                    NotImplemented(); 
            end           
        end
        
        function obj = getIndexSafety(obj,ObjStatePartition)
            
            if ~(strcmp(obj.TypeOfVectorField,'TCL') || strcmp(obj.TypeOfVectorField,'ChainInt'))
                error('This function can only be used for safety specifications. Please check if this is your case.')
            end
            
            % Retuns the indices of the safe set associated with the partition
            
            tempStatePartition = ObjStatePartition.Partition;
            
            
            switch obj.TypeOfVectorField
                case 'TCL'
                    NumberOfPoints = size(tempStatePartition.X,1);
                case 'ChainInt'
                    NumberOfPoints = size(tempStatePartition.grid_x,1);
                otherwise
                    NotImplemented();
            end
            
            
            % Information about the safe and reach sets
            SafeSet = obj.param.SafeSet;
            
            % Initializing the variables that will store the indices
            safeIndex = false(NumberOfPoints,1);
            
            % Iterating over states
            for i = 1:NumberOfPoints
                x = tempStatePartition.grid_x(i);  % current state
                
                % if current state belongs to safe set, then set the
                % corresponding index to true
                if min(SafeSet(1) <= x) && min(x <= SafeSet(2))
                    safeIndex(i) = true;
                end
            end
            
            % saving the results
            obj.IndexSafeSet = safeIndex;
            
        end
        
        function obj = BackwardIteration(obj,StatePartitionObj,InputPartition)
            
            OuterLoopInfo = obj.param.OuterLoopInfo;
            CurrentAmbiguity = obj.AmbiguityType;
            
            % Testing the value of N
            if isempty(obj.N)
                error('Please, initialize the field N before calling this function');
            elseif obj.N < 0 || ~isinteger(obj.N)
                error('N must be a positive integer (int8, int16, etc...)');
            end
            
            switch obj.TypeOfVectorField
                
                case 'TCL'
                    temp = PerformBackIteration(StatePartitionObj,InputPartition,obj.TypeOfVectorField,...
                        obj.AmbiguityType,obj.getIndexSafety(StatePartitionObj.getValues),...
                        @obj.iterateValueFunction,CurrentAmbiguity,OuterLoopInfo,obj.param);
                    
                case 'ChainInt'
                    
                    temp = PerformBackIteration(StatePartitionObj,InputPartition,obj.TypeOfVectorField,...
                        obj.AmbiguityType,obj.getIndexSafety(StatePartitionObj.getValues),...
                        @obj.iterateValueFunction,CurrentAmbiguity,OuterLoopInfo,obj.param);
                
                case 'Fishery'
                    
                    temp = PerformBackIteration(StatePartitionObj,InputPartition,obj.TypeOfVectorField,...
                        obj.AmbiguityType,obj.getIndexReachAvoid(StatePartitionObj.getValues),@obj.iterateValueFunction,...
                        CurrentAmbiguity,OuterLoopInfo,obj.param);
                    
                case 'CarPole'
                    
                    temp = PerformBackIteration(StatePartitionObj,InputPartition,obj.TypeOfVectorField,...
                        obj.AmbiguityType,obj.getIndexReachAvoid(StatePartitionObj.getValues),@obj.iterateValueFunction,...
                        CurrentAmbiguity,OuterLoopInfo,obj.param);
                    
                case 'CarPoleNL'
                    
                    temp = PerformBackIteration(StatePartitionObj,InputPartition,obj.TypeOfVectorField,...
                        obj.AmbiguityType,obj.getIndexReachAvoid(StatePartitionObj.getValues),@obj.iterateValueFunction,...
                        CurrentAmbiguity,OuterLoopInfo,obj.param);
                    
                otherwise
                    NotImplemented();
            end
            
            obj.ValueFunction = temp.ValueFunction;
            obj.OptInput = temp.OptInput;
            
            if isfield(temp,'ValueFunctionConservative')
                obj.ValueFunctionConservative = temp.ValueFunctionConservative;
                obj.OptInputConservative = temp.OptInputConservative;

                obj.ValueFunctionQP = temp.ValueFunctionQP;
                obj.OptInputQP = temp.OptInputQP;
            end
            
        end
        
        function out = iterateValueFunction(obj,CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel)
            
            StringCurrentState = createXandUString(CurrentState,[]);
            
            switch obj.TypeOfVectorField
                
                case 'TCL'
                    indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
                    
                    indexSafety = obj.IndexSafeSet;
                    
                    if ~indexSafety(indexCurrentState)
                        out = 0;
                    else
                        out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
                    end
                    
                case 'ChainInt'
                    
                    indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
                    
                    indexSafety = obj.IndexSafeSet;
                    
                    if ~indexSafety(indexCurrentState)
                        out = 0;
                    else
                        out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
                    end
                
                case 'Fishery'
                    
                    indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
                    
                    indexReachSet = obj.IndexSafeAndReachSet.reachIndex;
                    
                    if indexReachSet(indexCurrentState)
                        out = 1;
                    else
                        out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
                    end 
                    
                case 'CarPole'
                    
                    indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
                    
                    indexReachSet = obj.IndexSafeAndReachSet.reachIndex;
                    
                    if indexReachSet(indexCurrentState)
                        out = 1;
                    else
                        out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
                    end
                    
                case 'CarPoleNL'
                    
                    indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
                    
                    indexReachSet = obj.IndexSafeAndReachSet.reachIndex;
                    
                    if indexReachSet(indexCurrentState)
                        out = 1;
                    else
                        out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
                    end
                    
                otherwise
                    NotImplemented();
            end
            
            
        end
        
        function out = InnerOptimization(obj,x,u,NextValueFunc,StatePartitionObj,TypeOfKernel)
            
            % Returns the value function for a given state-action pair (x,u) using the
            % value function at the next iteration (NextValueFunc).
            
            TransitionProb = obj.param.TransitionProb; % vector with the transition probability matrix
            ListX = obj.param.ListX;
            
            if length(NextValueFunc) ~= length(ListX)
                error('The dimension of the thrid argument is incorrect.') % outputs an error if L is inconsistent with the size of the input Z
            end
            
            
            StringXU = createXandUString(x,u);
            
            if ~TransitionProb.isKey({StringXU}) % outputs an error is there is no element in TransitionProb with the label (stringX,stringU)
                error('The input-action pair is not a member of the transition probability')
            end          
            
            ObjFunc = NextValueFunc;
            TransXU = TransitionProb.values({StringXU});
            
            switch obj.AmbiguityType
                case 'NoAmbiguity'
                    out = TransXU{1}'*ObjFunc; % value function at the current state-action pair
                case 'MomentAmbiguity'
                    TempPartition = StatePartitionObj.getValues.Partition;
                    
                    %Ambiguity parameters
                    rhoMu = obj.ParamAmbiguity.rhoMu;
                    rhoSigma = obj.ParamAmbiguity.rhoSigma;
                    
                    [SupportSet,mu,Sigma] = ComputeSupportSetMuSigma(TransXU{1},TempPartition.grid_x,'WithSigma',obj.TypeOfVectorField);
                    
                    OptPro = MomentBasedAmbiguity(ObjFunc,Sigma,mu,rhoSigma,rhoMu,SupportSet);
                    out = AnalyseResults(OptPro,obj.AmbiguityType,[],TypeOfKernel);
                case 'WassersteinAmbiguity'
                    CenterBall = TransXU{1};
                    
                    %Ambiguity parameters
                    ep = obj.ParamAmbiguity.ep;
                    
                    OptPro = WassersteinAmbiguity(ObjFunc,ep,CenterBall);
                    
                    out = AnalyseResults(OptPro,obj.AmbiguityType,[],TypeOfKernel);
                case 'KernelAmbiguity'
                    CenterBall = TransXU{1};
                    gridX = StatePartitionObj.getValues.Partition.grid_x;
                    
                    OptPro = KernelBasedAmbiguity(ObjFunc,obj.ParamAmbiguity.ep,CenterBall,@GaussianKernel,obj.ParamAmbiguity.gamma,gridX);
                    
                    % Ambiguity parameterrs
                    OptPro.m = obj.ParamAmbiguity.m;
                    
                    OptPro.CurrentState = x;
                    OptPro.Input = u;  OptPro.param = obj.param;
                    
                    out = AnalyseResults(OptPro,obj.AmbiguityType,StatePartitionObj,TypeOfKernel);
                case 'KLdivAmbiguity'
                    
                    %Ambiguity parameters
                    ep = obj.ParamAmbiguity.ep;
                    
                    CenterBall = TransXU{1};
                    OptPro = KLdivAmbiguity(ObjFunc,ep,CenterBall);
                    out = AnalyseResults(OptPro,obj.AmbiguityType,[],TypeOfKernel);
                otherwise
                    error('This type of ambiguity has not been implemented. Please change the field AmbiguityType to a valid type.');
            end
            
            
        end
    end
end

function [safeIndex,reachIndex] = GetReachAvoid(dimensionStateSpace,StatePartition,param)

switch dimensionStateSpace
    case '3D'
        
        Nx1 = size(StatePartition.X1,2);
        Nx2 = size(StatePartition.X1,1);
        Nx3 = size(StatePartition.X1,3);
        
        % Information about the safe and reach sets
        SafeSet = param.SafeSet;
        ReachSet = param.ReachSet;
        
        % Initializing the variables that will store the indices
        safeIndex = false(Nx1*Nx2*Nx3,1);
        reachIndex = false(Nx1*Nx2*Nx3,1);
        
        % Iterating over states
        for i1 = 1:Nx1
            for i2=1:Nx2
                for i3=1:Nx3
                    x = [StatePartition.X1(i2,i1,i3);StatePartition.X2(i2,i1,i3);StatePartition.X3(i2,i1,i3)];  % current state
                    index = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],1,[]);
                    
                    % if current state belongs to safe set, then set the
                    % corresponding index to true
                    if min(SafeSet(:,1) <= x) && min(x <= SafeSet(:,2))
                        safeIndex(index) = true;
                    else
                        safeIndex(index) = false;
                    end
                    
                    % do the same with the reach set
                    if min(ReachSet(:,1) <= x) && min(x <= ReachSet(:,2))
                        reachIndex(index) = true;
                    else
                        reachIndex(index) = false;
                    end
                end
            end
        end
        
    case '4D'
        
        Nx1 = size(StatePartition.X1,1);
        Nx2 = size(StatePartition.X1,2);
        Nx3 = size(StatePartition.X1,3);
        Nx4 = size(StatePartition.X1,4);
        
        % Information about the safe and reach sets
        SafeSet = param.SafeSet;
        ReachSet = param.ReachSet;
        
        % Initializing the variables that will store the indices
        safeIndex = false(Nx1*Nx2*Nx3*Nx4,1);
        reachIndex = false(Nx1*Nx2*Nx3*Nx4,1);
        
        % Iterating over states
        for i1 = 1:Nx1
            for i2=1:Nx2
                for i3=1:Nx3
                    for i4=1:Nx4
                        x = [StatePartition.X1(i1,i2,i3,i4);StatePartition.X2(i1,i2,i3,i4);...
                            StatePartition.X3(i1,i2,i3,i4);StatePartition.X4(i1,i2,i3,i4)];  % current state
                        index = RemainingIterations(4,[[i1;i2;i3;i4],[Nx1;Nx2;Nx3;Nx4]],1,[]);
                        
                        % if current state belongs to safe set, then set the
                        % corresponding index to true
                        if min(SafeSet(:,1) <= x) && min(x <= SafeSet(:,2))
                            safeIndex(index) = true;
                        else
                            safeIndex(index) = false;
                        end
                        
                        % do the same with the reach set
                        if min(ReachSet(:,1) <= x) && min(x <= ReachSet(:,2))
                            reachIndex(index) = true;
                        else
                            reachIndex(index) = false;
                        end
                    end
                end
            end
        end
        
    otherwise
        NotImplemented();
end


end

function [SupportSet,mu,Sigma] = ComputeSupportSetMuSigma(Prob,Partition,type,TypeOfVectorField)

switch TypeOfVectorField
    case 'TCL'
        SupportSet = Partition';
    case 'Fishery'
        SupportSet = Partition';
    case 'ChainInt'
        SupportSet = Partition';
    case 'CarPole'
        SupportSet = Partition';
    case 'CarPoleNL'
        SupportSet = Partition';
    otherwise
        NotImplemented();
end
        
mu = SupportSet*Prob;
        
M = 1000;

switch TypeOfVectorField
    case 'TCL'
        if strcmp(type,'WithSigma')
            samples = discretesample(Prob,M);
            sumSigma = zeros(1,1);
            
            for i =1:M
                x = SupportSet(:,samples(i));
                temp = x - mu;
                sumSigma = temp*temp' + sumSigma;
            end
            
            Sigma = sumSigma/M;
        else
            Sigma = [];
        end
        
    case 'ChainInt'
        if strcmp(type,'WithSigma')
            samples = discretesample(Prob,M);
            sumSigma = zeros(2,2);
            
            for i =1:M
                x = SupportSet(:,samples(i));
                temp = x - mu;
                sumSigma = temp*temp' + sumSigma;
            end
            
            Sigma = sumSigma/M;
        else
            Sigma = [];
        end
        
    case 'Fishery'
        if strcmp(type,'WithSigma')
            samples = discretesample(Prob,M);
            sumSigma = zeros(3,3);
            
            for i =1:M
                x = SupportSet(:,samples(i));
                temp = x - mu;
                sumSigma = temp*temp' + sumSigma;
            end
            
            Sigma = sumSigma/M;
        else
            Sigma = [];
        end
        
    case 'CarPole'
        if strcmp(type,'WithSigma')
            samples = discretesample(Prob,M);
            sumSigma = zeros(4,4);
            
            for i =1:M
                x = SupportSet(:,samples(i));
                temp = x - mu;
                sumSigma = temp*temp' + sumSigma;
            end
            
            Sigma = sumSigma/M;
        else
            Sigma = [];
        end
        
    case 'CarPoleNL'
        if strcmp(type,'WithSigma')
            samples = discretesample(Prob,M);
            sumSigma = zeros(4,4);
            
            for i =1:M
                x = SupportSet(:,samples(i));
                temp = x - mu;
                sumSigma = temp*temp' + sumSigma;
            end
            
            Sigma = sumSigma/M;
        else
            Sigma = [];
        end
        
    otherwise
        NotImplemented();
end
        
end

function value = AnalyseResults(ObjAmbiguity,type,StatePartitionObj,TypeOfKernel)

if strcmp(type,'KernelAmbiguity')
    
    switch TypeOfKernel
        case 'Conservative'
            ObjAmbiguity = ObjAmbiguity.SolveOptimizationConservative(StatePartitionObj);
            value = ObjAmbiguity.OptRes.opt_value;
        case 'Matrix'
            ObjAmbiguity = ObjAmbiguity.SolveOptimization(StatePartitionObj);
            value = ObjAmbiguity.OptRes.opt_value;
        case 'QP'
            ObjAmbiguity = ObjAmbiguity.SolveOptimizationQP(StatePartitionObj);
            value = ObjAmbiguity.OptRes.opt_value;
        otherwise
            error('This type of Kernel has not been implemented.')
    end    
else
    ObjAmbiguity = ObjAmbiguity.SolveOptimization;
    if (ObjAmbiguity.OptRes.SolverStatus.problem == 0) || (ObjAmbiguity.OptRes.SolverStatus.problem == 4) || (ObjAmbiguity.OptRes.SolverStatus.problem == -1)
        value = ObjAmbiguity.OptRes.opt_value;
    else
        error('There is a problem when solving the optimization problem')
    end
end

end

function out = PerformBackIteration(StatePartitionObj,InputPartition,TypeOfVectorField,...
        AmbiguityType,IndexSafeAndOrReachSet,iterateValueFunction,...
        CurrentAmbiguity,OuterLoopInfo,param)

TimeHorizon = double(param.N);
NumberInputs = size(InputPartition,1); % Number of all possible combination of control inputs


switch TypeOfVectorField
    
     case 'TCL'
        
        NumberOfPoints = param.NumberOfPartitions; % number of points of the value function
        NumberInputs = size(InputPartition,1); % Number of all possible combination of control inputs
        
        
        if isempty(IndexSafeAndOrReachSet)
            error('This function cannot run if IndexSafeSet is empty. Please try using the method getIndexSafety');
        end
        
        
        ValueFunction = zeros(NumberOfPoints,TimeHorizon + 1);
        OptInput = zeros(NumberOfPoints,TimeHorizon + 1);
        ValueFunction(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        
        if strcmp(AmbiguityType,'KernelAmbiguity')
            ValueFunctionConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionConservative(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set

            ValueFunctionQP = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputQP = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionQP(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        end
        
    case 'ChainInt'
        
        NumberOfPoints = param.NumberOfPartitions(1)*param.NumberOfPartitions(2); % number of points of the value function
        NumberInputs = size(InputPartition,1); % Number of all possible combination of control inputs
        
        
        if isempty(IndexSafeAndOrReachSet)
            error('This function cannot run if IndexSafeSet is empty. Please try using the method getIndexSafety');
        end
        
        ValueFunction = zeros(NumberOfPoints,TimeHorizon + 1);
        OptInput = zeros(NumberOfPoints,TimeHorizon + 1);
        ValueFunction(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        
        if strcmp(AmbiguityType,'KernelAmbiguity')
            ValueFunctionConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionConservative(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set

            ValueFunctionQP = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputQP = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionQP(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        end
    
    case 'Fishery'
        NumberOfPoints = param.NumberOfPartitions(1)*param.NumberOfPartitions(2)*param.NumberOfPartitions(3); % number of points of the value function
        
        if isempty(IndexSafeAndOrReachSet)
            error('This function cannot run if IndexSafeAndReachSet is empty. Please try using the method getIndexReachAvoid');
        end
        
        
        ValueFunction = zeros(NumberOfPoints,TimeHorizon + 1);
        OptInput = zeros(NumberOfPoints,TimeHorizon + 1);
        ValueFunction(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        
        if strcmp(AmbiguityType,'KernelAmbiguity')
            ValueFunctionConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionConservative(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set

            ValueFunctionQP = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputQP = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionQP(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        end
        
    case 'CarPole'
        
        NumberOfPoints = param.NumberOfPartitions(1)*param.NumberOfPartitions(2)*param.NumberOfPartitions(3)*param.NumberOfPartitions(4); % number of points of the value function
        
        if isempty(IndexSafeAndOrReachSet)
            error('This function cannot run if IndexSafeAndReachSet is empty. Please try using the method getIndexReachAvoid');
        end
        
        
        ValueFunction = zeros(NumberOfPoints,TimeHorizon + 1);
        OptInput = zeros(NumberOfPoints,TimeHorizon + 1);
        ValueFunction(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set

        if strcmp(AmbiguityType,'KernelAmbiguity')
            ValueFunctionConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionConservative(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set

            ValueFunctionQP = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputQP = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionQP(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        end

    case 'CarPoleNL'
        
        NumberOfPoints = param.NumberOfPartitions(1)*param.NumberOfPartitions(2)*param.NumberOfPartitions(3)*param.NumberOfPartitions(4); % number of points of the value function
        
        if isempty(IndexSafeAndOrReachSet)
            error('This function cannot run if IndexSafeAndReachSet is empty. Please try using the method getIndexReachAvoid');
        end
        
        
        ValueFunction = zeros(NumberOfPoints,TimeHorizon + 1);
        OptInput = zeros(NumberOfPoints,TimeHorizon + 1);
        ValueFunction(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        
        if strcmp(AmbiguityType,'KernelAmbiguity')
            ValueFunctionConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputConservative = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionConservative(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set

            ValueFunctionQP = zeros(NumberOfPoints,TimeHorizon + 1);
            OptInputQP = zeros(NumberOfPoints,TimeHorizon + 1);
            ValueFunctionQP(IndexSafeAndOrReachSet.IndexSafeSet,TimeHorizon+1) = 1; % initializing the value function on the reach set
        end
        
    otherwise
        NotImplemented();
end



Grid_x = StatePartitionObj.getValues.Partition.grid_x;


% Creating a progress bar of the value function computation
total_iterations = TimeHorizon*(NumberOfPoints)*NumberInputs;
PrintInnerLoop(total_iterations,0,0,CurrentAmbiguity,OuterLoopInfo);

for i = TimeHorizon:-1:1
    NextValueFunc = ValueFunction(:,i+1); % saving in a temporary variable the value function of the next step
    ValueFunctionTemp = zeros(NumberOfPoints,1);
    OptInputTemp = zeros(NumberOfPoints,1);
    
    % Since we now have two value function for
    % the Kernel ambiguity we need to modify
    % this
    if strcmp(AmbiguityType,'KernelAmbiguity')
        NextValueFuncConservative = ValueFunctionConservative(:,i+1); % saving in a temporary variable the value function of the next step
        ValueFunctionTempConservative = zeros(NumberOfPoints,1);
        OptInputTempConservative = zeros(NumberOfPoints,1);

        NextValueFuncQP = ValueFunctionQP(:,i+1); % saving in a temporary variable the value function of the next step
        ValueFunctionTempQP = zeros(NumberOfPoints,1);
        OptInputTempQP = zeros(NumberOfPoints,1);
    end
    
    
    
    tic;

    for j = 1:NumberOfPoints % iterates over the number of points
        
        tempAmbiguityType = AmbiguityType; % This is necessary due to the parfor loop below

        x = Grid_x(j,:)'; % getting the current state to be updated
        
        % Since we now have two value function for
        % the Kernel ambiguity we need to modify
        % this
        if strcmp(tempAmbiguityType,'KernelAmbiguity')
            tempValueFuncMatrix = zeros(NumberInputs,1);
            tempValueFuncConservative = zeros(NumberInputs,1);
            tempValueFuncQP = zeros(NumberInputs,1);
        else
            tempValueFunc = zeros(NumberInputs,1);
        end
        
        for uCounter = 1:NumberInputs
            u = InputPartition(uCounter,:)'; % iterating over the number of inputs
            
            % Since we now have three value functions for
            % the Kernel ambiguity we need to modify
            % this
            if strcmp(tempAmbiguityType,'KernelAmbiguity')
                tempValueFuncConservative(uCounter)  = iterateValueFunction(x,u,NextValueFuncConservative,StatePartitionObj,'Conservative');	 % getting the new value for the value function
                tempValueFuncMatrix(uCounter) = iterateValueFunction(x,u,NextValueFunc,StatePartitionObj,'Matrix');	 % getting the new value for the value function
                tempValueFuncQP(uCounter) = iterateValueFunction(x,u,NextValueFuncQP,StatePartitionObj,'QP');	 % getting the new value for the value function
            else
                tempValueFunc(uCounter) = iterateValueFunction(x,u,NextValueFunc,StatePartitionObj,[]);	 % getting the new value for the value function
            end
        end
        
        if strcmp(tempAmbiguityType,'KernelAmbiguity')
            [ValueFunctionTemp(j),OptInputTemp(j)] = max(tempValueFuncMatrix); % storing the optimal value function and policy
            [ValueFunctionTempConservative(j),OptInputTempConservative(j)] = max(tempValueFuncConservative); % storing the optimal value function and policy
            [ValueFunctionTempQP(j),OptInputTempQP(j)] = max(tempValueFuncQP); % storing the optimal value function and policy
        else
            [ValueFunctionTemp(j),OptInputTemp(j)] = max(tempValueFunc); % storing the optimal value function and policy
        end

    end

    Lastime = toc;
    
    ValueFunction(:,i) = ValueFunctionTemp;
    OptInput(:,i) = OptInputTemp;
    
    if strcmp(AmbiguityType,'KernelAmbiguity')
        ValueFunctionConservative(:,i) = ValueFunctionTempConservative;
        OptInputConservative(:,i) = OptInputTempConservative;

        ValueFunctionQP(:,i) = ValueFunctionTempQP;
        OptInputQP(:,i) = OptInputTempQP;
    end
    % Printing on the screen
    tempInt = double(TimeHorizon-i+1);
    iterates = RemainingIterations(1,[tempInt,TimeHorizon],NumberOfPoints*NumberInputs,[]); % This is the number of iterations completes so far. The name of the matlab function may be misleading
    PrintInnerLoop(total_iterations,iterates,Lastime,CurrentAmbiguity,OuterLoopInfo);
end

out.ValueFunction = ValueFunction;
out.OptInput = OptInput;

if strcmp(AmbiguityType,'KernelAmbiguity')
    out.ValueFunctionConservative = ValueFunctionConservative;
    out.OptInputConservative = OptInputConservative;

    out.ValueFunctionQP = ValueFunctionQP;
    out.OptInputQP = OptInputQP;
end



end




