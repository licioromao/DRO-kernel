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
        function obj = KernelBasedAmbiguity(ObjFunc,ep,CenterBall)
            obj = obj@DistanceBasedAmbiguity(ObjFunc,ep,CenterBall); % Calling the connstructor of the parent class
            obj.gamma = [];
            obj.m = [];
            obj.param = [];
            obj.CurrentState = [];
            obj.Input = [];
        end
                       
        function obj = SolveOptimization(obj,ObjPartition)
            if isempty(obj.gamma) || isempty(obj.m) || isempty(obj.CurrentState) || isempty(obj.Input)
                error('You must first initialize the gamma and m fiels of the object, as well as CurrentState and Input using the SetValue function');
            end
            
            GaussianKernel = @(x,y) exp(-obj.gamma*norm(x-y)^2); % This is too specific for the Gaussian Kernel
            
         
            SumKernel = 0;
            alpha = obj.q;
            tempPartition = ObjPartition.getValues.Partition;
                       
            L = size(obj.param.List,1);
            for i=1:L
                for j=1:L
                    if i ~= L && j ~= L
                        SumKernel = SumKernel + alpha(i)*alpha(j)*GaussianKernel(tempPartition.grid_x(i,:)',tempPartition.grid_x(j,:)');
                    else
                        if i == L && j ~=L
                            SumKernel = SumKernel + alpha(i)*alpha(j)*GaussianKernel(-ones(3,1),tempPartition.grid_x(j,:)');
                        end
                        
                        if i == L && j == L
                            SumKernel = SumKernel + alpha(i)*alpha(j)*GaussianKernel(-ones(3,1),-ones(3,1));
                        end
                        
                        if i ~= L && j ==L
                            SumKernel = SumKernel + alpha(i)*alpha(j)*GaussianKernel(tempPartition.grid_x(i,:)',-ones(3,1));
                        end                        
                    end
                end
            end
                      
            obj.OptRes.opt_value = max(0,alpha'*obj.c - (obj.epsilon)*sqrt(SumKernel));
            
            %fprintf('Test prob: %.2f\n\n',(tempXSum.*alpha)'*obj.c - (obj.epsilon)*sqrt(SumKernel))
            
            %             if obj.OptRes >= 1
            %                 error('Error on the computation of the value function. It cannot be larger than one!')
            %             end
        end
    end
end


