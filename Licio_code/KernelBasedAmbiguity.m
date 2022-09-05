classdef KernelBasedAmbiguity < DistanceBasedAmbiguity
    % This class implements the computation of the Kernel-based ambiguity
    % set problem
    
    % This is very specific of the Gaussian Kernel. Can be generalized later.    
    properties
        param
        gamma
        m
        CurrentState
        Input
    end
       
    methods
        function obj = KernelBasedAmbiguity(arg1,arg2,arg3)
            obj = obj@DistanceBasedAmbiguity(arg1,arg2,arg3); % Calling the connstructor of the parent class
            obj.gamma = [];
            obj.m = [];
            obj.param = [];
            obj.CurrentState = [];
            obj.Input = [];
        end
                       
        function obj = SolveOptimization(obj)
            if isempty(obj.gamma) || isempty(obj.m) || isempty(obj.CurrentState) || isempty(obj.Input)
                error('You must first initialize the gamma and m fiels of the object, as well as CurrentState and Input using the SetValue function');
            end
            
            GaussianKernel = @(x,y) exp(-obj.gamma*norm(x-y)^2); % This is too specific for the Gaussian Kernel
            
            TempGrid = StatePartition(obj.param.NumberOfPartitions,obj.param.SafeSet);
            SamplesX = zeros(3,obj.m);
            KernelMatrix = zeros(obj.m,obj.m);
            alpha = obj.q;
            tempXSum = zeros(length(obj.param.List),1);
            
            
            for i=1:obj.m               
                Noise = generateNoise(obj.param);
                TempIndex = TempGrid.computeElementPartition(obj.CurrentState,obj.Input,Noise,obj.param);
                NextStatePartition = TempGrid.getValues.grid_x(TempIndex.elementPartition);
                CenterNextState = TempIndex.nextState;
                
                if isempty(CenterNextState)
                    SamplesX(:,i) = -ones(3,1);
                else
                    SamplesX(:,i) = CenterNextState(1:3);
                end
                
                key = BuildString(NextStatePartition);
                tempXSum = tempXSum + FindIndex(key,obj.param.List,alpha);
            end
            
            for i=1:obj.m
                for j=i:obj.m
                    KernelMatrix(i,j) = GaussianKernel(SamplesX(:,i),SamplesX(:,j)); 
                end
            end
            
            KernelMatrix = (KernelMatrix + KernelMatrix')/2;
                   
            
            obj.OptRes.opt_value = min(1,tempXSum'*obj.c + (1/obj.m)^2*sum(sum(KernelMatrix))*obj.epsilon);
            
%             if obj.OptRes >= 1
%                 error('Error on the computation of the value function. It cannot be larger than one!')
%             end
        end
    end
end

function out = BuildString(key)

if isempty(key)
    out = 'NaN';
else
    out = sprintf('(%.2f,%.2f,%.2f)',key);
end

end


function out = FindIndex(key,List,alpha)


L = length(List);
goal = false;
i = 1;

out = zeros(L,1);

while ~goal && i <= L
    if strcmp(key,List{i})
        out(i) = alpha(i);
        goal = true;
    end
    i = i + 1;
end

if i == L+2
    error('Key not found. There should be an error in the code.')
end

end

