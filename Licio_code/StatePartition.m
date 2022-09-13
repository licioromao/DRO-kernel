classdef StatePartition
       
    properties (Access = private)
        Partition
        TypeOfVectorField
    end
    
    methods
        function obj = StatePartition(NumberOfPartitions,SafeSet,TypeOfVectorField)
            obj.TypeOfVectorField = TypeOfVectorField;
            obj.Partition = generatePartition(NumberOfPartitions,SafeSet,obj.TypeOfVectorField);
        end
        
        function out = getValues(obj)
            out.Partition = obj.Partition;
            out.TypeOfVectorField = obj.TypeOfVectorField;
        end
        
        function [out,outX] = createList(obj,InputPartition)
            
            % This function creates a vector of array that contains the labels of the
            % states of the generated discrete model
            
            % Inputs: StatePartition -- this is the output of the function
            % generatePartition as defined above
            %
            % Ouput: out -- labels for the discrete states
            
            switch obj.TypeOfVectorField
                case 'Fishery'         
                    Nx1 = size(obj.Partition.X1,2);
                    Nx2 = size(obj.Partition.X2,1);
                    Nx3 = size(obj.Partition.X3,3);
                    Nu = size(InputPartition,1);
                    
                    out = cell(Nx1*Nx2*Nx3*Nu + 1,1);
                    outX = cell(Nx1*Nx2*Nx3 + 1,1);
                    
                    for i1 = 1:Nx1
                        for i2 = 1:Nx2
                            for i3 = 1:Nx3
                                x = [obj.Partition.X1(i2,i1,i3),obj.Partition.X2(i2,i1,i3),obj.Partition.X3(i2,i1,i3)];
                                for j=1:Nu
                                    u = InputPartition(j,:);
                                    
                                    tempIndex1 = RemainingIterations(4,[[i1;i2;i3;j],[Nx1;Nx2;Nx3;Nu]],1,[]);
                                    
                                    out{tempIndex1} = createXandUString(x,u);
                                end
                                
                                tempIndex2 = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],1,[]);
                                outX{tempIndex2} = createXandUString(x,[]);
                            end
                        end
                    end
                    
                    out{end} = 'NaN'; % The last state called 'NaN' represents a fictious state of the discrete model containing unsafe states
                    outX{end} = 'NaN'; % The last state called 'NaN' represents a fictious state of the discrete model containing unsafe states
                case 'TCL'
                    
                    N = size(obj.Partition.X,1);
                    Nu = 2;
                    
                    out = cell(Nu*N+1,1);
                    outX = cell(N + 1,1);
                    
                    for i = 1:N
                        x = obj.Partition.X(i);
                        for j=1:2
                            
                            tempIndex1 = RemainingIterations(2,[[i;j],[N;2]],1,[]);
                            
                            out{tempIndex1} = createXandUString(x,j-1);
                        end
                        outX{i} = createXandUString(x,[]);
                    end
                    
                    out{end} = 'NaN'; % The last state called 'NaN' represents a fictious state of the discrete model containing unsafe states
                    outX{end} = 'NaN'; % The last state called 'NaN' represents a fictious state of the discrete model containing unsafe states
                otherwise
                    NotImplemented();
            end
            
        end
        
        function out = getSizePartition(obj)
            
            % This function returns the granularity of the partition in each dimension.
            % Its input is the output of the function generatePartition.
            switch obj.TypeOfVectorField
                case 'Fishery'
                    h1 = obj.Partition.X1(1,2,1) - obj.Partition.X1(1,1,1);
                    h2 = obj.Partition.X2(2,1,1) - obj.Partition.X2(1,1,1);
                    h3 = obj.Partition.X3(1,1,2) - obj.Partition.X3(1,1,1);
                    
                    out = [h1;h2;h3];
                case 'TCL'
                    out = obj.Partition.X(2) - obj.Partition.X(1);
                otherwise
                    NotImplemented();
            end
        end
        
        function out = getElementPartition(obj,x)
            
            % This function returns the element of the partition that x belongs to. If
            % the input x is outside the safe set, this function returns an empty
            % array.
            switch obj.TypeOfVectorField
                case 'Fishery'
                    index = [find(obj.Partition.X1(1,:,1) <= x(1),1,'last');
                        find(obj.Partition.X2(:,1,1) <= x(2),1,'last');
                        find(obj.Partition.X3(1,1,:) <= x(3),1,'last')]; % Check the element of the partition
                    
                    if length(index) ~= 3 % If there is one dimensional for which index is empty, set out to empty array
                        index = [];
                        out.index = [];
                        out.x = [];
                    end
                    
                    if ~isempty(index) % Otherwise, the function outputs a structure with index and xHat field
                        xHat = [obj.Partition.X1(1,index(1),1);obj.Partition.X2(index(2),1,1);obj.Partition.X3(1,1,index(3));x(4)];
                        out.index = index;
                        out.x = xHat;
                    end
                case 'TCL'
                    index = find(obj.Partition.X <= x,1,'last');
                    
                    if isempty(index)
                        out.index = [];
                        out.x = [];
                    else
                        xHat = obj.Partition.X(index);
                        out.index = index;
                        out.x = xHat;
                    end
                    
                otherwise
                    NotImplemented();
            end
            
        end
        
        function out = getCenterPartition(obj,x)
            
            % This function receives as input a point x in the lattice and returns the
            % center of the partition.
            
            switch obj.TypeOfVectorField
                case 'Fishery'
                    out = zeros(4,1);
                    
                    if isempty(x)
                        out = [];
                    else
                        out(1:3) = x(1:3) + obj.getSizePartition/2;
                        out(4) = x(4);
                    end
                case 'TCL'
                                     
                    if isempty(x)
                        out = [];
                    else
                        out = x + obj.getSizePartition/2;
                    end
                    
                otherwise
                    NotImplemented();
            end
        end
        
        function out = computeElementPartition(obj,currentState,input,noise,param)
            
            % Given currentState, input and noise vectors, this function uses the
            % vector_field class to progate the dynamics and returns the element and
            % center of the partition corresponding to the next state.
            
            out.nextState = [];
            out.elementPartition = [];
            
            switch obj.TypeOfVectorField
                case 'Fishery'
                    
                    VF = VectorFieldFishery(currentState,input,param); % creating a vector field object
                    VF.Noise = noise; % setting the noise parameter in the vector field class
                    nextState = VF.IterateDynamics; % generating next state
                    
                case 'TCL'
                    
                    VF = VectorFieldTCL(currentState,input,param);
                    VF.Noise = noise;
                    nextState = VF.IterateDynamics;
                    
                otherwise
                    NotImplemented();
            end
            
            temp = obj.getElementPartition(nextState);
            out.elementPartition = temp.index;
            out.nextState = obj.getCenterPartition(temp.x);
        end
        
    end
end

