classdef StatePartition
       
    properties (Access = private)
        Partition
    end
    
    methods
        function obj = StatePartition(NumberOfPartitions,SafeSet)
            obj.Partition = generatePartition(NumberOfPartitions,SafeSet);
        end
        
        function out = getValues(obj)
            out = obj.Partition;
        end
        
        function out = createList(obj)
            
            % This function creates a vector of array that contains the labels of the
            % states of the generated discrete model
            
            % Inputs: StatePartition -- this is the output of the function
            % generatePartition as defined above
            %
            % Ouput: out -- labels for the discrete states
            
            Nx1 = size(obj.Partition.X1,2);
            Nx2 = size(obj.Partition.X2,1);
            Nx3 = size(obj.Partition.X3,3);
            
            out = {};
            
            for i1 = 1:Nx1
                for i2 = 1:Nx2
                    for i3 = 1:Nx3
                        x = [obj.Partition.X1(i2,i1,i3),obj.Partition.X2(i2,i1,i3),obj.Partition.X3(i2,i1,i3)];
                        out = [out;sprintf('(%.2f,%.2f,%.2f)',x(1),x(2),x(3))];
                    end
                end
            end
            
            out = [out;'NaN']; % The last state called 'NaN' represents a fictious state of the discrete model containing unsafe states
            
        end
        
        function out = getSizePartition(obj)
            
            % This function returns the granularity of the partition in each dimension.
            % Its input is the output of the function generatePartition.
            
            h1 = obj.Partition.X1(1,2,1) - obj.Partition.X1(1,1,1);
            h2 = obj.Partition.X2(2,1,1) - obj.Partition.X2(1,1,1);
            h3 = obj.Partition.X3(1,1,2) - obj.Partition.X3(1,1,1);
            
            out = [h1;h2;h3];
        end
        
        function out = getElementPartition(obj,x)
            
            % This function returns the element of the partition that x belongs to. If
            % the input x is outside the safe set, this function returns an empty
            % array.
            
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
            
        end
        
        function out = getCenterPartition(obj,x)
            
            % This function receives as input a point x in the lattice and returns the
            % center of the partition.
            
            out = zeros(4,1);
            
            if isempty(x)
                out = [];
            else
                out(1:3) = x(1:3) + obj.getSizePartition/2;
                out(4) = x(4);
            end
            
        end
        
        function out = computeElementPartition(obj,currentState,input,noise,param)
            
            % Given currentState, input and noise vectors, this function uses the
            % vector_field class to progate the dynamics and returns the element and
            % center of the partition corresponding to the next state.
            
            out.nextState = [];
            out.elementPartition = [];
            
            VF = VectorFieldFishery(currentState,input,param); % creating a vector field object
            VF.Noise = noise; % setting the noise parameter in the vector field class
            nextState = VF.IterateDynamics; % generating next state
            
            
            temp = obj.getElementPartition(nextState);
            out.elementPartition = temp.index;
            out.nextState = obj.getCenterPartition(temp.x);
        end
        
    end
end

