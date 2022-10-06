classdef StatePartition
       
    properties (Access = private)
        Partition
        TypeOfVectorField
    end
    
    methods
        function obj = StatePartition(NumberOfPartitions,SafeSet,TypeOfVectorField)
            obj.TypeOfVectorField = TypeOfVectorField;
            switch obj.TypeOfVectorField
                case 'TCL'
                    obj.Partition = generatePartition(NumberOfPartitions,SafeSet,'1D');
                case 'ChainInt'
                    obj.Partition = generatePartition(NumberOfPartitions,SafeSet,'2D');
                case 'Fishery'
                    obj.Partition = generatePartition(NumberOfPartitions,SafeSet,'3D');
                case 'CarPole'
                    obj.Partition = generatePartition(NumberOfPartitions,SafeSet,'4D');
                case 'CarPoleNL'
                    obj.Partition = generatePartition(NumberOfPartitions,SafeSet,'4D');
                otherwise
                    NotImplemented();
            end
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
                case 'TCL'
                    
                    [out,outX] = ObjCreateList(obj.Partition,InputPartition,'1D');
                    
                case 'ChainInt'
                    
                    [out,outX] = ObjCreateList(obj.Partition,InputPartition,'2D');
                                        
                case 'Fishery'         
                    
                    [out,outX] = ObjCreateList(obj.Partition,InputPartition,'3D');
                    
                 case 'CarPole'         
                    
                    [out,outX] = ObjCreateList(obj.Partition,InputPartition,'4D');
                    
                case 'CarPoleNL'
                     
                    [out,outX] = ObjCreateList(obj.Partition,InputPartition,'4D');
                     
                otherwise
                    NotImplemented();
            end
            
        end
        
        function out = getSizePartition(obj)
            
            % This function returns the granularity of the partition in each dimension.
            % Its input is the output of the function generatePartition.
            switch obj.TypeOfVectorField
                case 'TCL'
                    out = obj.Partition.X(2) - obj.Partition.X(1);
                    
                case 'ChainInt'
                    h1 = obj.Partition.X1(1,2) - obj.Partition.X1(1,1);
                    h2 = obj.Partition.X2(2,1) - obj.Partition.X2(1,1);
                    
                    out = [h1;h2];
                    
                case 'Fishery'
                    h1 = obj.Partition.X1(1,2,1) - obj.Partition.X1(1,1,1);
                    h2 = obj.Partition.X2(2,1,1) - obj.Partition.X2(1,1,1);
                    h3 = obj.Partition.X3(1,1,2) - obj.Partition.X3(1,1,1);
                    
                    out = [h1;h2;h3];
                    
                case 'CarPole'
                    
                    h1 = obj.Partition.X1(2,1,1,1) - obj.Partition.X1(1,1,1,1);
                    h2 = obj.Partition.X2(1,2,1,1) - obj.Partition.X2(1,1,1,1);
                    h3 = obj.Partition.X3(1,1,2,1) - obj.Partition.X3(1,1,1,1);
                    h4 = obj.Partition.X4(1,1,1,2) - obj.Partition.X4(1,1,1,1);
                    
                    out = [h1;h2;h3;h4];
                    
                case 'CarPoleNL'
                    
                    h1 = obj.Partition.X1(2,1,1,1) - obj.Partition.X1(1,1,1,1);
                    h2 = obj.Partition.X2(1,2,1,1) - obj.Partition.X2(1,1,1,1);
                    h3 = obj.Partition.X3(1,1,2,1) - obj.Partition.X3(1,1,1,1);
                    h4 = obj.Partition.X4(1,1,1,2) - obj.Partition.X4(1,1,1,1);
                    
                    out = [h1;h2;h3;h4];
                
                otherwise
                    NotImplemented();
            end
        end
        
        function out = getElementPartition(obj,x)
            
            % This function returns the element of the partition that x belongs to. If
            % the input x is outside the safe set, this function returns an empty
            % array.
            switch obj.TypeOfVectorField
                
                case 'TCL'
                    index = find(obj.Partition.X <= x,1,'last');
                    
                    if isempty(index)
                        out.index = 1;
                        out.x = obj.Partition.X(1);
                    else
                        xHat = obj.Partition.X(index);
                        out.index = index;
                        out.x = xHat;
                    end
                case 'ChainInt'
                    
                    index = [find(obj.Partition.X1(1,:) <= x(1),1,'last');
                        find(obj.Partition.X2(:,1) <= x(2),1,'last');];
                    
                    if length(index) ~= 2
                        out.index = [1;1];
                        out.x = obj.Partition.grid_x(1,:)';
                    else
                        xHat = [obj.Partition.X1(1,index(1));obj.Partition.X2(index(2),1)];
                        out.index = index;
                        out.x = xHat;
                    end
                
                case 'Fishery'
                    index = [find(obj.Partition.X1(1,:,1) <= x(1),1,'last');
                        find(obj.Partition.X2(:,1,1) <= x(2),1,'last');
                        find(obj.Partition.X3(1,1,:) <= x(3),1,'last')]; % Check the element of the partition
                    
                    if length(index) ~= 3 % If there is one dimensional for which index is empty, set out to empty array
                        out.index = [1;1;1];
                        out.x = obj.Partition.grid_x(1,:)';
                        out.x = [out.x;1];
                    else
                        xHat = [obj.Partition.X1(1,index(1),1);obj.Partition.X2(index(2),1,1);obj.Partition.X3(1,1,index(3));x(4)];
                        out.index = index;
                        out.x = xHat;
                    end
                    
                case 'CarPole'
                    index = [find(obj.Partition.X1(:,1,1,1) <= x(1),1,'last');
                        find(obj.Partition.X2(1,:,1,1) <= x(2),1,'last');
                        find(obj.Partition.X3(1,1,:,1) <= x(3),1,'last'); 
                        find(obj.Partition.X3(1,1,1,:) <= x(4),1,'last')]; % Check the element of the partition
                    
                    if length(index) ~= 4 % If there is one dimensional for which index is empty, set out to empty array
                        out.index = [1;1;1;1];
                        out.x = obj.Partition.grid_x(1,:)';
                        out.x = out.x;
                    else
                        xHat = [obj.Partition.X1(index(1),1,1,1);obj.Partition.X2(1,index(2),1,1);obj.Partition.X3(1,1,index(3),1);obj.Partition.X4(1,1,1,index(4))];
                        out.index = index;
                        out.x = xHat;
                    end
                    
                case 'CarPoleNL'
                    
                    index = [find(obj.Partition.X1(:,1,1,1) <= x(1),1,'last');
                        find(obj.Partition.X2(1,:,1,1) <= x(2),1,'last');
                        find(obj.Partition.X3(1,1,:,1) <= x(3),1,'last');
                        find(obj.Partition.X3(1,1,1,:) <= x(4),1,'last')]; % Check the element of the partition
                    
                    if length(index) ~= 4 % If there is one dimensional for which index is empty, set out to empty array
                        out.index = [1;1;1;1];
                        out.x = obj.Partition.grid_x(1,:)';
                        out.x = out.x;
                    else
                        xHat = [obj.Partition.X1(index(1),1,1,1);obj.Partition.X2(1,index(2),1,1);obj.Partition.X3(1,1,index(3),1);obj.Partition.X4(1,1,1,index(4))];
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
                case 'TCL'  
                    out = x + obj.getSizePartition/2; 
                case 'ChainInt'
                    out = x + obj.getSizePartition/2;
                case 'Fishery'
                    out = zeros(4,1);
                    out(1:3) = x(1:3) + obj.getSizePartition/2;
                    out(4) = x(4);
                case 'CarPole'
                    out = x + obj.getSizePartition/2;
                case 'CarPoleNL'
                    out = x + obj.getSizePartition/2;
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
                case 'TCL'
                    
                    VF = VectorFieldTCL(currentState,input,param);
                    VF.Noise = noise;
                    nextState = VF.IterateDynamics;
                    
                case 'ChainInt'
                    
                    VF = VectorFieldChainInt(currentState,input,param);
                    VF.Noise = noise;
                    nextState = VF.IterateDynamics;
                    
                case 'Fishery'
                    
                    VF = VectorFieldFishery(currentState,input,param); % creating a vector field object
                    VF.Noise = noise; % setting the noise parameter in the vector field class
                    nextState = VF.IterateDynamics; % generating next state
                    
                case 'CarPole'
                    
                    VF = VectorFieldCarPole(currentState,input,param);
                    VF.Noise = noise;
                    nextState = VF.IterateDynamics;
                    
                case 'CarPoleNL'
                    
                    VF = VectorFieldCarPoleNL(currentState,input,param);
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


function [List,ListX] = ObjCreateList(StatePartition,InputPartition,dimensionState)

switch dimensionState
    case '1D'
        N = size(StatePartition.X,1);
        Nu = size(InputPartition,1);
        
        List = cell(Nu*N,1);
        ListX = cell(N,1);
        
        for i = 1:N
            x = StatePartition.X(i);
            for j=1:Nu
                
                tempIndex1 = RemainingIterations(2,[[i;j],[N;2]],1,[]);
                
                List{tempIndex1} = createXandUString(x,InputPartition(j));
            end
            ListX{i} = createXandUString(x,[]);
        end
        
    case '2D'
        
        Nx1 = size(c.X1,2);
        Nx2 = size(StatePartition.X2,1);
        
        Nu = size(InputPartition,1);
        
        List = cell(Nx1*Nx2*Nu,1);
        ListX = cell(Nx1*Nx2,1);
        
        for i1 = 1:Nx1
            for i2 = 1:Nx2
                
                x = [obj.Partition.X1(i2,i1),obj.Partition.X2(i2,i1)];
                for j=1:Nu
                    u = InputPartition(j,:);
                    
                    tempIndex1 = RemainingIterations(3,[[i1;i2;j],[Nx1;Nx2;Nu]],1,[]);
                    
                    List{tempIndex1} = createXandUString(x,u);
                end
                
                tempIndex2 = RemainingIterations(2,[[i1;i2],[Nx1;Nx2]],1,[]);
                ListX{tempIndex2} = createXandUString(x,[]);
                
            end
        end
        
    case '3D'
        
        Nx1 = size(StatePartition.X1,2);
        Nx2 = size(StatePartition.X2,1);
        Nx3 = size(StatePartition.X3,3);
        Nu = size(InputPartition,1);
        
        List = cell(Nx1*Nx2*Nx3*Nu,1);
        ListX = cell(Nx1*Nx2*Nx3,1);
        
        for i1 = 1:Nx1
            for i2 = 1:Nx2
                for i3 = 1:Nx3
                    x = [StatePartition.X1(i2,i1,i3),StatePartition.X2(i2,i1,i3),StatePartition.X3(i2,i1,i3)];
                    for j=1:Nu
                        u = InputPartition(j,:);
                        
                        tempIndex1 = RemainingIterations(4,[[i1;i2;i3;j],[Nx1;Nx2;Nx3;Nu]],1,[]);
                        
                        List{tempIndex1} = createXandUString(x,u);
                    end
                    
                    tempIndex2 = RemainingIterations(3,[[i1;i2;i3],[Nx1;Nx2;Nx3]],1,[]);
                    ListX{tempIndex2} = createXandUString(x,[]);
                end
            end
        end
        
    case '4D'
        
        Nx1 = size(StatePartition.X1,1);
        Nx2 = size(StatePartition.X1,2);
        Nx3 = size(StatePartition.X1,3);
        Nx4 = size(StatePartition.X1,4);
        
        Nu = size(InputPartition,1);
        
        List = cell(Nx1*Nx2*Nx3*Nx4*Nu,1);
        ListX = cell(Nx1*Nx2*Nx3*Nx4,1);
        
        for i1 = 1:Nx1
            for i2 = 1:Nx2
                for i3 = 1:Nx3
                    for i4 = 1:Nx4
                        x = [StatePartition.X1(i1,i2,i3,i4),StatePartition.X2(i1,i2,i3,i4),StatePartition.X3(i1,i2,i3,i4),StatePartition.X4(i1,i2,i3,i4)];
                        for j=1:Nu
                            u = InputPartition(j,:);
                            
                            tempIndex1 = RemainingIterations(5,[[i1;i2;i3;i4;j],[Nx1;Nx2;Nx3;Nx4;Nu]],1,[]);
                            
                            List{tempIndex1} = createXandUString(x,u);
                        end
                        
                        tempIndex2 = RemainingIterations(4,[[i1;i2;i3;i4],[Nx1;Nx2;Nx3;Nx4]],1,[]);
                        ListX{tempIndex2} = createXandUString(x,[]);
                    end
                end
            end
        end
        
    otherwise
        NotImplemented();     
end

end



