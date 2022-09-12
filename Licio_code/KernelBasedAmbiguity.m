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
                       
        function obj = SolveOptimization(obj,ObjPartition)
            if isempty(obj.gamma) || isempty(obj.m) || isempty(obj.CurrentState) || isempty(obj.Input)
                error('You must first initialize the gamma and m fiels of the object, as well as CurrentState and Input using the SetValue function');
            end
            
            GaussianKernel = @(x,y) exp(-obj.gamma*norm(x-y)^2); % This is too specific for the Gaussian Kernel
            
          %  SamplesX = zeros(3,obj.m);
            SumKernel = 0;
            alpha = obj.q;
            tempXSum = zeros(length(obj.param.List),1);
            tempPartition = ObjPartition.getValues.Partition;
            
%             for i=1:obj.m
%                 Noise = generateNoise(obj.param);
%                 TempIndex = ObjPartition.computeElementPartition(obj.CurrentState,obj.Input,Noise,obj.param);
%                                 
%                 if isempty(TempIndex.elementPartition)
%                     SamplesX(:,i) = -ones(3,1);
%                     NextStatePartition = [];
%                 else
%                     SamplesX(:,i) = TempIndex.nextState(1:3);
%                     NextStatePartition = [tempPartition.X1(1,TempIndex.elementPartition(1),1),tempPartition.X2(TempIndex.elementPartition(2),1,1),tempPartition.X3(1,1,TempIndex.elementPartition(3))];
%                 end
%                 
%                 key = BuildString(NextStatePartition);
%                 tempXSum = tempXSum + FindIndex(key,obj.param.List);
%             end
%             
%             tempXSum = tempXSum/obj.m;
%             
%             fprintf('Norm transition:%.2f\n',norm(tempXSum/obj.m - alpha))
%             if norm(tempXSum/obj.m - alpha) > 0.4
%                 fprintf('\n\n\nx=(%.2f,%.2f,%.2f) u=(%.2f,%.2f)\n\n\n',obj.CurrentState,obj.Input);
%                 
%                 Noise = generateNoise(obj.param);
%                 TempIndex = ObjPartition.computeElementPartition(obj.CurrentState,obj.Input,Noise,obj.param);
%                 NextStatePartition = ObjPartition.getValues.grid_x(TempIndex.elementPartition);
%                 
%                 error('Stop Here!')
%             end
            
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

function out = BuildString(index)

if isempty(index)
    out = 'NaN';
else
    out = sprintf('(%.2f,%.2f,%.2f)',index);
end

end


function out = FindIndex(key,List)


L = length(List);
goal = false;
i = 1;

out = zeros(L,1);

while ~goal && i <= L
    if strcmp(key,List{i})
        out(i) = 1;
        goal = true;
    end
    i = i + 1;
end

if i == L+1 && ~goal
    error('Key not found. There should be an error in the code.')
end

end
% 
% function ProbMeas = FindTransitionProb(x,u,TransitionMatrix)
% 
% key =false;
% L = length(TransitionMatrix);
% StringX = sprintf('(%.2f,%.2f,%.2f)',x);
% StringU = sprintf('(%.2f,%.2f)',u);
% i= 1;
% 
% while ~key && i <= L
%     if strcmp(TransitionMatrix{i}.State,StringX) && strcmp(TransitionMatrix{i}.Action,StringU)
%         key = true;
%         ProbMeas = TransitionMatrix{i}.ProbMeasure;
%     end
%     i = i + 1;
% end
% 
% end

