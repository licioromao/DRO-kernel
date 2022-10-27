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
        KernelFunc
        CholKernelMatrix
    end
       
    methods
        function obj = KernelBasedAmbiguity(ObjFunc,ep,CenterBall,FuncKernel,gamma,grid_x)
            obj = obj@DistanceBasedAmbiguity(ObjFunc,ep,CenterBall); % Calling the connstructor of the parent class
            obj.gamma = gamma;
            obj.m = [];
            obj.param = [];
            obj.CurrentState = [];
            obj.Input = [];
            obj.KernelFunc = FuncKernel;
            obj.CholKernelMatrix = FuncKernel(grid_x',grid_x',gamma,[],'MatrixFac');
        end
                       
        function obj = SolveOptimizationConservative(obj,ObjPartition)
            if isempty(obj.gamma) || isempty(obj.m) || isempty(obj.CurrentState) || isempty(obj.Input)
                error('You must first initialize the gamma and m fiels of the object, as well as CurrentState and Input using the SetValue function');
            end

            alpha = obj.q;
            grid_x = ObjPartition.getValues.Partition.grid_x;

            tempKernel = obj.KernelFunc(grid_x',grid_x',obj.gamma,alpha,'Normal');

            obj.OptRes.opt_value = max(0,alpha'*obj.c - (obj.epsilon)*tempKernel-2*obj.epsilon);
        end

        function obj = SolveOptimization(obj,ObjPartition)
            if isempty(obj.gamma) || isempty(obj.m) || isempty(obj.CurrentState) || isempty(obj.Input)
                error('You must first initialize the gamma and m fiels of the object, as well as CurrentState and Input using the SetValue function');
            end


            alpha = obj.q;
            grid_x = ObjPartition.getValues.Partition.grid_x;

            beta = obj.CholKernelMatrix\(obj.CholKernelMatrix'\obj.c);

            tempKernel = obj.KernelFunc(grid_x',grid_x',obj.gamma,beta,'Normal');

            obj.OptRes.opt_value = max(0,alpha'*obj.c - (obj.epsilon)*tempKernel);

        end

        function obj = SolveOptimizationQP(obj,ObjPartition)

            if isempty(obj.gamma) || isempty(obj.m) || isempty(obj.CurrentState) || isempty(obj.Input)
                error('You must first initialize the gamma and m fiels of the object, as well as CurrentState and Input using the SetValue function');
            end

            grid_x = ObjPartition.getValues.Partition.grid_x;
            M = length(grid_x); %number of sample points 

            alpha = sdpvar(M,1);
            Constraints = [sum(alpha) == 1, alpha >= 0]; % alpha is a probability vector

            tempKernel = obj.CholKernelMatrix; tempKernel = tempKernel'*tempKernel;  % Kernel matrix
            
            Constraints = [Constraints,alpha'*tempKernel*alpha - 2/M*alpha'*tempKernel*ones(M,1) + (1/M^2)*ones(M,1)'*tempKernel*ones(M,1) <= obj.epsilon^2];
            
            options = sdpsettings('solver','mosek','verbose',0); % options to solve the optimisation problem
            objective = -obj.c'*alpha; % objective function

            sol = optimize(Constraints,objective,options);

            % Error if unfeasible
            if sol.problem == 1
                error('Unfeasible optimization problem')
            end
            
            % Saving the results
            
            obj.OptRes = [];
            obj.OptRes.SolverStatus = sol; % information about opt problem
            %obj.OptRes.dualVar = dual(Constraints);
            
            obj.OptRes.Nvar = length(depends(Constraints)); % number of optimisation variable (this needs to be doubled-checked)
            
            obj.OptRes.p = value(alpha);
            %fprintf('In this formulation, we cannot recover the value of the optimal distribution \n');
            
            obj.OptRes.opt_value = value(-objective);

        end
        
    end
end


