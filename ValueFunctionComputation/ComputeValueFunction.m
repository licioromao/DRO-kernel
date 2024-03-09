classdef ComputeValueFunction
    % NEED TO ADD A DETAILED DESCRIPTION HERE
    
    properties
        
        value_function % This value stores the optimal value function
        opt_input % This value stores the optimal action

        time_horizon % horizon of the reach-avoid property
        time % this will store the time to go through a full value function computation
    end
    
    methods
        function obj = ComputeValueFunction(number_of_points,time_horizon)

            obj.value_function = zeros(number_of_points,time_horizon+1); % initializing the variable to store the value function
            obj.time_horizon = time_horizon;
            obj.time = [];

        end
    end
    
    methods (Abstract)
        backward_iteration(obj,state_partition,input_partition)

        iterate_value_function(obj,current_state,current_input,...
                                next_value_func,state_partition,type_of_kernel)

        inner_optimisation(obj,x,u,next_value_func,state_partition,...
                            type_of_kernel)
    end
   
end        
%         function obj = getIndexReachAvoid(obj,ObjStatePartition)
%             
%             if ~(strcmp(obj.type_vector_field,'Fishery') || strcmp(obj.type_vector_field,'CarPole') || strcmp(obj.type_vector_field,'CarPoleNL') )
%                 error('This function can only be used for reach avoid specifications. Please check if this is your case.')
%             end
%             
%             % Retuns the indices of the safe and reach sets associated with the partition
%             
%             tempStatePartition = ObjStatePartition.Partition;
%             
%             switch obj.type_vector_field
%                 case 'Fishery'
%                     [obj.index_safe_reach_set.safeIndex,obj.index_safe_reach_set.reachIndex] = GetReachAvoid('3D',tempStatePartition,obj.param);
%                 case 'CarPole'
%                     [obj.index_safe_reach_set.safeIndex,obj.index_safe_reach_set.reachIndex] = GetReachAvoid('4D',tempStatePartition,obj.param);
%                 case 'CarPoleNL'
%                     [obj.index_safe_reach_set.safeIndex,obj.index_safe_reach_set.reachIndex] = GetReachAvoid('4D',tempStatePartition,obj.param);
%                 otherwise
%                     NotImplemented(); 
%             end           
%         end
        
        
        
%         function out = iterateValueFunction(obj,CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel)
%             
%             StringCurrentState = createXandUString(CurrentState,[]);
%             
%             switch obj.type_vector_field
%                 
%                 case 'TCL'
%                     indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
%                     
%                     indexSafety = obj.index_safe_set;
%                     
%                     if ~indexSafety(indexCurrentState)
%                         out = 0;
%                     else
%                         out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
%                     end
%                     
%                 case 'ChainInt'
%                     
%                     indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
%                     
%                     indexSafety = obj.index_safe_set;
%                     
%                     if ~indexSafety(indexCurrentState)
%                         out = 0;
%                     else
%                         out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
%                     end
%                 
%                 case 'Fishery'
%                     
%                     indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
%                     
%                     indexReachSet = obj.index_safe_reach_set.reachIndex;
%                     
%                     if indexReachSet(indexCurrentState)
%                         out = 1;
%                     else
%                         out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
%                     end 
%                     
%                 case 'CarPole'
%                     
%                     indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
%                     
%                     indexReachSet = obj.index_safe_reach_set.reachIndex;
%                     
%                     if indexReachSet(indexCurrentState)
%                         out = 1;
%                     else
%                         out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
%                     end
%                     
%                 case 'CarPoleNL'
%                     
%                     indexCurrentState = strcmp(obj.param.ListX,StringCurrentState);
%                     
%                     indexReachSet = obj.index_safe_reach_set.reachIndex;
%                     
%                     if indexReachSet(indexCurrentState)
%                         out = 1;
%                     else
%                         out = obj.InnerOptimization(CurrentState,Input,NextValueFunc,StatePartitionObj,TypeOfKernel);
%                     end
%                     
%                 otherwise
%                     NotImplemented();
%             end
%             
%             
%         end
        
%         function out = InnerOptimization(obj,x,u,NextValueFunc,StatePartitionObj,TypeOfKernel)

%                 case 'KLdivAmbiguity'
%                     
%                     %Ambiguity parameters
%                     ep = obj.param_ambiguity.ep;
%                     
%                     CenterBall = TransXU{1};
%                     OptPro = KLdivAmbiguity(ObjFunc,ep,CenterBall);
%                     out = AnalyseResults(OptPro,obj.ambiguity_type,[],TypeOfKernel);
%                 otherwise
%                     error('This type of ambiguity has not been implemented. Please change the field AmbiguityType to a valid type.');
%             end
%             
%             
%         end

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





