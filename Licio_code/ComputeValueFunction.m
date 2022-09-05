classdef ComputeValueFunction
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ValueFunction % This value stores the optimal value function
        OptInput % This value stores the optimal action
        AmbiguityType % Type of ambiguity set
        param % structure with all the parameters of the problem
        IndexSafeAndReachSet % Indices of the safe and reach set
        N % horizon of the reach-avoid property
    end
    
    methods
        function obj = ComputeValueFunction(ArgParam,Type)
            NumberOfPoints = ArgParam.NumberOfPartitions(1)*ArgParam.NumberOfPartitions(2)*ArgParam.NumberOfPartitions(3); % number of points of the value function
            obj.ValueFunction = zeros(NumberOfPoints + 1,ArgParam.N+1); % initializing the variable to store the value function
            obj.OptInput = zeros(NumberOfPoints + 1,ArgParam.N+1); % initializing the variable to store the optimal input
            obj.AmbiguityType = Type;
            obj.param = ArgParam;
            obj.IndexSafeAndReachSet = [];
            obj.N = ArgParam.N;
        end
        
        function obj = getIndexReachAvoid(obj,StatePartition)
            
            % Retuns the indices of the safe and reach sets associated with the partition
            
            Nx1 = size(StatePartition.X1,2);
            Nx2 = size(StatePartition.X1,1);
            Nx3 = size(StatePartition.X1,3);
            
            % Information about the safe and reach sets
            SafeSet = obj.param.SafeSet;
            ReachSet = obj.param.ReachSet;
            
            % Initializing the variables that will store the indices
            safeIndex = false(Nx1*Nx2*Nx3 +1,1);
            reachIndex = false(Nx1*Nx2*Nx3 + 1,1);
            
            % Iterating over states
            for i1 = 1:Nx1
                for i2=1:Nx2
                    for i3=1:Nx3
                        x = [StatePartition.X1(i2,i1,i3);StatePartition.X2(i2,i1,i3);StatePartition.X3(i2,i1,i3)];  % current state
                        index = (i1-1)*Nx2*Nx3 + (i2-1)*Nx3 + i3;
                        
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
            
            % The last state correspond to the unsafe state of the MDP, so we set to
            % false
            safeIndex(end) = false;
            reachIndex(end) = false;
            
            
            % saving the results
            obj.IndexSafeAndReachSet.safeIndex = safeIndex;
            obj.IndexSafeAndReachSet.reachIndex = reachIndex;
            
        end
        
        function obj = BackwardIteration(obj,StatePartitionObj,InputPartition)
            
            % Testing the value of N
            if isempty(obj.N)
                error('Please, initialize the field N before calling this function');
            elseif obj.N < 0 || ~isinteger(obj.N)
                error('N must be a positive integer (int8, int16, etc...)');
            end
            
            NumberOfPoints = obj.param.NumberOfPartitions(1)*obj.param.NumberOfPartitions(2)*obj.param.NumberOfPartitions(3); % number of points of the value function
            NumberInputs = size(InputPartition,1); % Number of all possible combination of control inputs
            obj.ValueFunction(obj.IndexSafeAndReachSet.reachIndex,obj.N+1) = 1; % initializing the value function on the reach set
            Grid_x = StatePartitionObj.getValues.grid_x;
                    
            % Creating a progress bar of the value function computation
            hh = waitbar(0,'Initializing','Name','Computing Value Function...');
            total_iterations = double((obj.N)*NumberOfPoints*NumberInputs);
            fprintf('Total of iterations: %d \n',total_iterations);
            
            for i = obj.N:-1:1
                NextValueFunc = obj.ValueFunction(:,i+1); % saving in a temporary variable the value function of the next step
                ValueFunctionTemp = zeros(NumberOfPoints,1);
                OptInputTemp = zeros(NumberOfPoints,1);
                for j = 1:NumberOfPoints % iterates over the number of points                   
                    x = Grid_x(j,:)'; % getting the current state to be updated
                    tempValueFunc = zeros(NumberInputs,1);
                    for uCounter = 1:NumberInputs
                        u = InputPartition(uCounter,:)'; % iterating over the number of inputs
                        tempValueFunc(uCounter) = obj.iterateValueFunction(x,u,NextValueFunc,StatePartitionObj); % getting the new value for the value function
                    end
                    
                    % Update waitbar
                    iterates = RemainingIterations(2,[[obj.N-i+1;j],[obj.N;NumberOfPoints]],NumberInputs,hh);
                    perc_iterates = iterates/total_iterations;
                    waitbar(perc_iterates,hh,sprintf('%.5f completed',perc_iterates));
                    
                    [ValueFunctionTemp(j),OptInputTemp(j)] = max(tempValueFunc); % storing the optimal value function and policy
                end
                
                obj.ValueFunction(1:end-1,i) = ValueFunctionTemp;
                obj.OptInput(1:end -1,i) = OptInputTemp;
            end
            close(hh)
        end
        
        function out = iterateValueFunction(obj,CurrentState,Input,NextValueFunc,StatePartitionObj)
            
            indexCurrentState = strcmp(obj.param.List,sprintf('(%.2f,%.2f,%.2f)',CurrentState(1),CurrentState(2),CurrentState(3)));
            
            indexReachSet = obj.IndexSafeAndReachSet.reachIndex;
            
            if indexReachSet(indexCurrentState)
                out = 1;
            else
                out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj);
            end
            
        end
        
        function out = InnerOptimization(obj,x,u,NextValueFunc,StatePartitionObj)
            
            % Returns the value function for a given state-action pair (x,u) using the
            % value function at the next iteration (NextValueFunc). 
            
            TransitionProb = obj.param.TransitionProb; % vector with the transition probability matrix
            L = length(TransitionProb); % number of states in the transition probability
            
            if length(NextValueFunc) ~= length(TransitionProb{1}.ProbMeasure)
                error('The dimension of the thrid argument is incorrect.') % outputs an error if L is inconsistent with the size of the input Z
            end
            
            stringX = sprintf('(%.2f,%.2f,%.2f)',x(1),x(2),x(3)); % creating the string of the current state
            stringU = sprintf('(%.2f,%.2f)',u(1),u(2)); % creating the string of the current action
            
            key = false;
            i = 1;
            TempPartition = StatePartition(obj.param.NumberOfPartitions,obj.param.SafeSet);
                   
            % Searching for the correct transition probability
            while i <= L && ~key
                if strcmp(TransitionProb{i}.State,stringX) && strcmp(TransitionProb{i}.Action,stringU) % if stringX and stringU matches with the information in the transition prob matrix
                    key = true;
                    ObjFunc = NextValueFunc;
                    switch obj.AmbiguityType
                        case 'NoAmbiguity'
                            out = TransitionProb{i}.ProbMeasure'*ObjFunc; % value function at the current state-action pair
                        case 'MomentAmbiguity'
                            [SupportSet,mu,Sigma] = ComputeSupportSetMuSigma(TransitionProb{i}.ProbMeasure,TempPartition.getValues.grid_x,'WithSigma');
                            rhoMu = 2; rhoSigma = 2;
                            OptPro = MomentBasedAmbiguity(ObjFunc,Sigma,mu,rhoMu,rhoSigma,SupportSet);  
                            out = AnalyseResults(OptPro,obj.AmbiguityType,[]);                        
                        case 'WassersteinAmbiguity'
                            ep = 0.1; CenterBall = TransitionProb{i}.ProbMeasure;
                            OptPro = WassersteinAmbiguity(ep,ObjFunc,CenterBall);
                            out = AnalyseResults(OptPro,obj.AmbiguityType,[]);
                        case 'KernelAmbiguity'
                            ep = 0.1; CenterBall = TransitionProb{i}.ProbMeasure;
                            OptPro = KernelBasedAmbiguity(ep,ObjFunc,CenterBall);
                            OptPro.gamma = 5;  OptPro.CurrentState = x;
                            OptPro.Input = u; OptPro.m = 50; OptPro.param = obj.param;
                            out = AnalyseResults(OptPro,obj.AmbiguityType,StatePartitionObj);
                        case 'KLdivAmbiguity'
                            ep = 0.01; CenterBall = TransitionProb{i}.ProbMeasure;
                            OptPro = KLdivAmbiguity(ep,ObjFunc,CenterBall);
                            out = AnalyseResults(OptPro,obj.AmbiguityType,[]);
                        otherwise
                            error('This type of ambiguity has not been implemented. Please change the field AmbiguityType to a valid type.');
                    end
                end
                i = i+1;
            end
            
            if ~key % outputs an error is there is no element in TransitionProb with the label (stringX,stringU)
                error('The input-action pair is not a member of the transition probability')
            end
            
        end


    end
end


function [SupportSet,mu,Sigma] = ComputeSupportSetMuSigma(Prob,Partition,type)

SupportSet = Partition';
SupportSet = [SupportSet,[-1;-1;-1]];

mu = SupportSet*Prob;

M = 1000;

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

end

function value = AnalyseResults(ObjAmbiguity,type,StatePartitionObj)

if strcmp(type,'KernelAmbiguity')
    ObjAmbiguity = ObjAmbiguity.SolveOptimization(StatePartitionObj);
    value = ObjAmbiguity.OptRes.opt_value;
else
    ObjAmbiguity = ObjAmbiguity.SolveOptimization;
    if (ObjAmbiguity.OptRes.SolverStatus.problem == 0) || (ObjAmbiguity.OptRes.SolverStatus.problem == 4)
        value = ObjAmbiguity.OptRes.opt_value;
    else
        error('There is a problem when solving the optimization problem')
    end
end

end



