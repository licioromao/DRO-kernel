classdef ComputeValueFunction
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ValueFunction % This value stores the optimal value function
        TypeOfVectorField % This value stores the type of Vector field
        OptInput % This value stores the optimal action
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
                case 'Fishery'
                    NumberOfPoints = ArgParam.NumberOfPartitions(1)*ArgParam.NumberOfPartitions(2)*ArgParam.NumberOfPartitions(3); % number of points of the value function
                case 'TCL'
                    NumberOfPoints = ArgParam.NumberOfPartitions; % number of points of the value function
                otherwise
                    NotImplemented();
            end
            
            
            obj.ValueFunction = zeros(NumberOfPoints,ArgParam.N+1); % initializing the variable to store the value function
            obj.OptInput = zeros(NumberOfPoints,ArgParam.N+1); % initializing the variable to store the optimal input
            obj.TypeOfVectorField = TypeOfVectorField;
            obj.AmbiguityType = StructAmbiguityTypes.Name;
            obj.ParamAmbiguity = StructAmbiguityTypes;
            obj.param = ArgParam;
            obj.IndexSafeAndReachSet = [];
            obj.N = ArgParam.N;
            obj.time = [];
        end
        
        function obj = getIndexReachAvoid(obj,ObjStatePartition)
            
            if ~strcmp(obj.TypeOfVectorField,'Fishery')
                error('This function can only be used for reach avoid specifications. Please check if this is your case.')
            end
            
            % Retuns the indices of the safe and reach sets associated with the partition
            
            tempStatePartition = ObjStatePartition.Partition;
            
            
            Nx1 = size(tempStatePartition.X1,2);
            Nx2 = size(tempStatePartition.X1,1);
            Nx3 = size(tempStatePartition.X1,3);
            
            % Information about the safe and reach sets
            SafeSet = obj.param.SafeSet;
            ReachSet = obj.param.ReachSet;
            
            % Initializing the variables that will store the indices
            safeIndex = false(Nx1*Nx2*Nx3,1);
            reachIndex = false(Nx1*Nx2*Nx3,1);
            
            % Iterating over states
            for i1 = 1:Nx1
                for i2=1:Nx2
                    for i3=1:Nx3
                        x = [tempStatePartition.X1(i2,i1,i3);tempStatePartition.X2(i2,i1,i3);tempStatePartition.X3(i2,i1,i3)];  % current state
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
            
            % saving the results
            obj.IndexSafeAndReachSet.safeIndex = safeIndex;
            obj.IndexSafeAndReachSet.reachIndex = reachIndex;
            
        end
        
        function obj = getIndexSafety(obj,ObjStatePartition)
            
            if ~strcmp(obj.TypeOfVectorField,'TCL')
                error('This function can only be used for safety specifications. Please check if this is your case.')
            end
            
            % Retuns the indices of the safe set associated with the partition
            
            tempStatePartition = ObjStatePartition.Partition;
            
            
            NumberOfPoints = size(tempStatePartition.X,1);
            
            % Information about the safe and reach sets
            SafeSet = obj.param.SafeSet;
            
            % Initializing the variables that will store the indices
            safeIndex = false(NumberOfPoints,1);
            
            % Iterating over states
            for i = 1:NumberOfPoints
                x = tempStatePartition.X(i);  % current state
                
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
                case 'Fishery'
                    
                    TimeHorizon = double(obj.N);
                    
                    NumberOfPoints = obj.param.NumberOfPartitions(1)*obj.param.NumberOfPartitions(2)*obj.param.NumberOfPartitions(3); % number of points of the value function
                    NumberInputs = size(InputPartition,1); % Number of all possible combination of control inputs
                    
                    if isempty(obj.IndexSafeAndReachSet)
                        error('This function cannot run if IndexSafeAndReachSet is empty. Please try using the method getIndexReachAvoid');
                    end
                    
                    obj.ValueFunction(obj.IndexSafeAndReachSet.reachIndex,obj.N+1) = 1; % initializing the value function on the reach set
                    Grid_x = StatePartitionObj.getValues.Partition.grid_x;
                    
%                    % Creating a progress bar of the value function computation
%                     hh = waitbar(0,'Initializing','Name','Computing Value Function...');
%                     total_iterations = TimeHorizon*(NumberOfPoints + 1)*NumberInputs;
%                     fprintf('Total of iterations: %d \n',total_iterations);
                    
                    % Creating a progress bar of the value function computation 
                    total_iterations = TimeHorizon*(NumberOfPoints)*NumberInputs;
                    
                    PrintInnerLoop(total_iterations,0,0,CurrentAmbiguity,OuterLoopInfo);
                    
                    for i = TimeHorizon:-1:1
                        NextValueFunc = obj.ValueFunction(:,i+1); % saving in a temporary variable the value function of the next step
                        ValueFunctionTemp = zeros(NumberOfPoints,1);
                        OptInputTemp = zeros(NumberOfPoints,1);
                        sumTime = 0;
                        sumIt = 0;
                        for j = 1:NumberOfPoints % iterates over the number of points
                            x = Grid_x(j,:)'; % getting the current state to be updated
                            tempValueFunc = zeros(NumberInputs,1);
                            
                            tic;
                            tempIterateFunc = @obj.iterateValueFunction;
                            parfor uCounter = 1:NumberInputs
                                u = InputPartition(uCounter,:)'; % iterating over the number of inputs
                                tempValueFunc(uCounter) = tempIterateFunc(x,u,NextValueFunc,StatePartitionObj); % getting the new value for the value function
                            end
                            Lastime = toc;
                            sumTime = sumTime + Lastime;
                            sumIt = sumIt + 1;
                            % Update waitbar
                            tempInt = TimeHorizon - i +1;
                            iterates = RemainingIterations(2,[[tempInt;j],[TimeHorizon;NumberOfPoints]],NumberInputs,[]); % This is the number of iterations completes so far. The name of the matlab function may be misleading
                            
                            if NumberOfPoints > 200
                                if mod(j,floor(0.01*NumberOfPoints)) == 0
                                    PrintInnerLoop(total_iterations,iterates,sumTime/sumIt,CurrentAmbiguity,OuterLoopInfo);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            else
                                if (mod(j,2) == 0)
                                    PrintInnerLoop(total_iterations,iterates,sumTime/sumIt,CurrentAmbiguity,OuterLoopInfo);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            end
                            
                            [ValueFunctionTemp(j),OptInputTemp(j)] = max(tempValueFunc); % storing the optimal value function and policy
                        end
                        
                        obj.ValueFunction(1:end,i) = ValueFunctionTemp;
                        obj.OptInput(1:end,i) = OptInputTemp;
                    end
                case 'TCL'
                    
                    TimeHorizon = double(obj.N);
                    
                    NumberOfPoints = obj.param.NumberOfPartitions; % number of points of the value function
                    NumberInputs = size(InputPartition,1); % Number of all possible combination of control inputs
                    
                    if isempty(obj.IndexSafeSet)
                        error('This function cannot run if IndexSafeAndReachSet is empty. Please try using the method getIndexSafety');
                    end                    
                    
                    obj.ValueFunction(obj.IndexSafeSet,obj.N+1) = 1; % initializing the value function on the reach set
                    Grid_x = StatePartitionObj.getValues.Partition.grid_x;
                    
                    % Creating a progress bar of the value function computation 
                    total_iterations = TimeHorizon*(NumberOfPoints)*NumberInputs;
                    PrintInnerLoop(total_iterations,0,0,CurrentAmbiguity,OuterLoopInfo);
                    
                    for i = TimeHorizon:-1:1
                        NextValueFunc = obj.ValueFunction(:,i+1); % saving in a temporary variable the value function of the next step
                        ValueFunctionTemp = zeros(NumberOfPoints,1);
                        OptInputTemp = zeros(NumberOfPoints,1);
                        
                        sumTime = 0;
                        sumIt = 0;
                        
                        for j = 1:NumberOfPoints % iterates over the number of points
                            x = Grid_x(j); % getting the current state to be updated
                            tempValueFunc = zeros(NumberInputs,1);
                            tic;
                            tempIterateFunc = @obj.iterateValueFunction;
                            parfor uCounter = 1:NumberInputs
                                u = InputPartition(uCounter); % iterating over the number of inputs
                                tempValueFunc(uCounter) = tempIterateFunc(x,u,NextValueFunc,StatePartitionObj);	 % getting the new value for the value function
                            end
                            Lastime = toc;
                            sumTime = sumTime + Lastime;
                            sumIt = sumIt + 1;
                            % Update waitbar
                            tempInt = double(TimeHorizon-i+1);
                            iterates = RemainingIterations(2,[[tempInt;j],[TimeHorizon;NumberOfPoints]],NumberInputs,[]); % This is the number of iterations completes so far. The name of the matlab function may be misleading
                            
                            if NumberOfPoints > 200
                                if mod(j,floor(0.01*NumberOfPoints)) == 0
                                    PrintInnerLoop(total_iterations,iterates,sumTime/sumIt,CurrentAmbiguity,OuterLoopInfo);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            else
                                if (mod(j,10) == 0)
                                    PrintInnerLoop(total_iterations,iterates,sumTime/sumIt,CurrentAmbiguity,OuterLoopInfo);
                                    sumTime = 0;
                                    sumIt = 0;
                                end
                            end
                            %waitbar(perc_iterates,hh,sprintf('%.5f completed. %.2f seconds to go',perc_iterates,SecToGo));
                            
                            [ValueFunctionTemp(j),OptInputTemp(j)] = max(tempValueFunc); % storing the optimal value function and policy
                        end
                        
                        obj.ValueFunction(1:end,i) = ValueFunctionTemp;
                        obj.OptInput(1:end,i) = OptInputTemp;
                    end
            end
        end
        
        function out = iterateValueFunction(obj,CurrentState,Input,NextValueFunc,StatePartitionObj)
            
            StringCurrentState = createXandUString(CurrentState,[]);
            
            switch obj.TypeOfVectorField
                case 'Fishery'
                    
                    indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
                    
                    indexReachSet = obj.IndexSafeAndReachSet.reachIndex;
                    
                    if indexReachSet(indexCurrentState)
                        out = 1;
                    else
                        out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj);
                    end
                case 'TCL'
                    indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
                    
                    indexSafety = obj.IndexSafeSet;
                    
                    if ~indexSafety(indexCurrentState)
                        out = 0;
                    else
                        out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj);
                    end
                    
                otherwise
                    NotImplemented();
            end
            
            
        end
        
        function out = InnerOptimization(obj,x,u,NextValueFunc,StatePartitionObj)
            
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
                    out = AnalyseResults(OptPro,obj.AmbiguityType,[]);
                case 'WassersteinAmbiguity'
                    CenterBall = TransXU{1};
                    
                    %Ambiguity parameters
                    ep = obj.ParamAmbiguity.ep;
                    
                    OptPro = WassersteinAmbiguity(ObjFunc,ep,CenterBall);
                    
                    out = AnalyseResults(OptPro,obj.AmbiguityType,[]);
                case 'KernelAmbiguity'
                    CenterBall = TransXU{1};
                    
                    OptPro = KernelBasedAmbiguity(ObjFunc,obj.ParamAmbiguity.ep,CenterBall,@GaussianKernel);
                    
                    % Ambiguity parameterrs
                    OptPro.gamma = obj.ParamAmbiguity.gamma;
                    OptPro.m = obj.ParamAmbiguity.m;
                    
                    OptPro.CurrentState = x;
                    OptPro.Input = u;  OptPro.param = obj.param;
                    
                    out = AnalyseResults(OptPro,obj.AmbiguityType,StatePartitionObj);
                case 'KLdivAmbiguity'
                    
                    %Ambiguity parameters
                    ep = obj.ParamAmbiguity.ep;
                    
                    CenterBall = TransXU{1};
                    OptPro = KLdivAmbiguity(ObjFunc,ep,CenterBall);
                    out = AnalyseResults(OptPro,obj.AmbiguityType,[]);
                otherwise
                    error('This type of ambiguity has not been implemented. Please change the field AmbiguityType to a valid type.');
            end
            
            
        end
    end
end


function [SupportSet,mu,Sigma] = ComputeSupportSetMuSigma(Prob,Partition,type,TypeOfVectorField)

switch TypeOfVectorField
    case 'TCL'
        SupportSet = Partition';
        
    case 'Fishery'
        SupportSet = Partition';
    otherwise
        NotImplemented();
end
        
mu = SupportSet*Prob;
        
M = 1000;

switch TypeOfVectorField
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
    otherwise
        NotImplemented();
end
        
end

function value = AnalyseResults(ObjAmbiguity,type,StatePartitionObj)

if strcmp(type,'KernelAmbiguity')
    ObjAmbiguity = ObjAmbiguity.SolveOptimization(StatePartitionObj);
    value = ObjAmbiguity.OptRes.opt_value;
else
    ObjAmbiguity = ObjAmbiguity.SolveOptimization;
    if (ObjAmbiguity.OptRes.SolverStatus.problem == 0) || (ObjAmbiguity.OptRes.SolverStatus.problem == 4) || (ObjAmbiguity.OptRes.SolverStatus.problem == -1)
        value = ObjAmbiguity.OptRes.opt_value;
    else
        error('There is a problem when solving the optimization problem')
    end
end

end



