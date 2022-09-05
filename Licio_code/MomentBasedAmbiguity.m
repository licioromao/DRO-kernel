classdef MomentBasedAmbiguity < Ambiguity
    % This class is a subclass of the abstract class ambiguity and
    % implements the abstract method SolveOptimization using a moment-based
    % formulation as in the Automatica paper by Insoon Yang "A dynamic game
    % approach to distributionally robust safety specifications for
    % stochastic systems". To define the ambiguity set using moment
    % information we follow the formulation of equation (5) in Insoon's
    % paper with slightly different notation.
    
    
    %  Comparison of the fields of the class with the notation of Insoon's
    %  paper.
    %
    %   * rhoSigma is the parameter defining the variance of the ambiguity
    %   set. It stands for the parameter c in equation (5) of Insoon's
    %   paper.
    %
    %   * rhoMu is the parameter defining the expectation radius in Insoon's
    %   paper. It stands for the parameter b in equation (5) of Insoon's
    %   paper.
    %
    %
       
    % Private fields of the class
    properties (Access = protected)
        mu  % "center" of the mean
        sigma % "center" of the variance
        rhoSigma % Radius of the "variance" characterizing the ambiguity set
        rhoMu % Radius of the "expectation" characterizing the ambiguity set
        supportSet % This is the discrete set from where the random variable takes value
    end
    
    methods
        % Constructor of the class
        function obj = MomentBasedAmbiguity(arg1,arg2,arg3,arg4,arg5,arg6)
          
            [ObjFunc,MeanCenter,VarianceCenter,SupportSet,rhoMu,rhoSigma] = verifyArg(arg1,arg2,arg3,arg4,arg5,arg6); % This function verifies the consistency of the parameters
            
            obj = obj@Ambiguity(ObjFunc); % Calls the constructor of the class ambiguity
            
            % Saving the values of the private field of the class
            obj.supportSet = SupportSet;
            obj.mu = MeanCenter;
            obj.sigma = VarianceCenter;
            
            obj.rhoMu = rhoMu;
            obj.rhoSigma = rhoSigma;       
            
        end
        
        function out = getValues(obj)
            
            % This function gets the values of the private fields and
            % returns them as a structure
            
            out.c = obj.c;
            out.mu = obj.mu;
            out.sigma = obj.sigma;
            out.rhoSigma = obj.rhoSigma;
            out.rhoMu = obj.rhoMu;
            out.supportSet = obj.supportSet;
        end
        
        function obj = setValues(obj,arg1,arg2,arg3,arg4,arg5,arg6)
            
            % This function modifies the private properties of the class.
            % It mimics the structure of its constructor 
            
            [ObjFunc,MeanCenter,VarianceCenter,SupportSet,RhoMu,RhoSigma] = verifyArg(arg1,arg2,arg3,arg4,arg5,arg6);
            
            % Updating the values of the private fields
            obj.c = ObjFunc;
            
            obj.supportSet = SupportSet;
            obj.mu = MeanCenter;
            obj.sigma = VarianceCenter;
            
            obj.rhoMu = RhoMu;
            obj.rhoSigma = RhoSigma; 
            
        end
        
        function obj = SolveOptimization(obj)
            
            % This function solves the optimization problem defined in the
            % description of the class Ambiguity. The formulation below is
            % based on the paper by Insoon Yang mentioned at the top of
            % this file
            
            m = size(obj.supportSet,2); % number of samples
            
            x = sdpvar(m,1); % optimization variable -> this will result in the optimal distribution
                
            Constraints = [];
            Constraints = [Constraints, x>=0, sum(x) == 1]; % setting the constraint for x to be a probability distribution
            % The next two constraints define the ambiguity set based on
            % moments
            Constraints = [Constraints, abs(obj.supportSet*x - obj.mu)  <= obj.rhoMu];
            temp = (obj.supportSet - obj.mu);
            Constraints = [Constraints, temp*diag(x)*temp'<= obj.rhoSigma*obj.sigma];
            
            options = sdpsettings('solver','mosek','verbose',0); % setting properties to solve the optimization variable
            
            sol = optimize(Constraints,obj.c'*x,options);
            
            % Error if unfeasible
            if sol.problem == 1
                error('Unfeasible optimization problem')
            end
            
            % Saving the results of the solution to the optimization
            % problems
            
            obj.OptRes = [];          
            obj.OptRes.SolverStatus = sol; % information about opt problem
            %obj.OptRes.dualVar = dual(Constraints);
            
            obj.OptRes.Nvar = length(depends(Constraints)); % number of optimisation variable (this needs to be doubled-checked)
            
            obj.OptRes.p = value(x);
            obj.OptRes.opt_value = value(obj.c'*x);
            
            obj.p = value(x);
        end
   
        function obj = SolveOptimizationDual(obj)
            
            % This function is similar to the methods SolveOptimization,
            % with the difference being that it finds the optimal solution
            % by solving the dual problem as described in Theorem 2 in
            % Insoon's paper.
            
            m = size(obj.supportSet,2); % Number of samples
            n = size(obj.supportSet,1); % dimension of the Euclidean space from which the random variable takes value
            
            % These are optimization variable according to the formulation
            % in the paper of Insoon
            lambda_bar = sdpvar(n,1); 
            bar_lambda = sdpvar(n,1);
            Lambda = sdpvar(n,n);
            nu = sdpvar(1);
                
            Constraints = [];
            
            Constraints = [Constraints,bar_lambda >= 0, lambda_bar >=0, Lambda >= 0];
            
            for i = 1:m
                Constraints = [Constraints, obj.c(i) + (obj.supportSet(:,i) - obj.mu)'*Lambda*(obj.supportSet(:,i) - obj.mu) + (obj.supportSet(:,i) - obj.mu)'*(lambda_bar - bar_lambda) + nu >= 0];
            end
            
            options = sdpsettings('solver','mosek','verbose',0);
            
            objective = bar_lambda'*ones(n,1)*obj.rhoMu + lambda_bar'*ones(n,1)*obj.rhoMu + nu + obj.rhoSigma*trace(Lambda*obj.sigma);
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
            
            obj.OptRes.p = [];
            %fprintf('In this formulation, we cannot recover the value of the optimal distribution \n');
            
            obj.OptRes.opt_value = value(-objective);
        end
    end
end

%% Auxiliary functions to verify the arguments
function [ObjFunc,MeanCenter,VarianceCenter,SupportSet,rhoMu,rhoSigma] = verifyArg(arg1,arg2,arg3,arg4,arg5,arg6)

% Setting the output to be empty arrays
MeanCenter = []; 
ObjFunc = arg1;
VarianceCenter = [];
SupportSet = [];
tol = 1e-5;

% Checking whether the first input is an array of numbers
if isvector(ObjFunc) && ~ischar(ObjFunc)
    % Converting ObjFunc into column vector if necessary
    if size(ObjFunc,1) == 1 
        ObjFunc = ObjFunc';
    end
else
    error('The first argument must be a row or column vector of real numbers')
end


% Checking whether there are two consecutive positive scalars as the arguments of the
% function. Then checks whether the remaing arguments has the structure
% (mu,Sigma,SupporSet) or (Sigma,mu,Support). An error appears if one of
% the conditions hold:
%
%   * The number of scalars in (arg2,...,arg6) is not equal to two
%   * Scalars are not consecutive
%   * There exists one negative scalar 
%   * There is no matrix in (arg2,...,arg6)
%   * There is no square or PSD or symmetric matrix in (arg2,...,arg6)
%   * Dimension of ObjFunc above is different than the number of columns in
%   the matrix that is not square, PSD or symmetric
%   * Dimension of the square and PSD matrix is not equal to the number of
%   rows of the SupportSet



if isscalar(arg2) && ~ischar(arg2) || isscalar(arg3) && ~ischar(arg3) || isscalar(arg4) && ~ischar(arg4) || isscalar(arg5) && ~ischar(arg5) ||  isscalar(arg6) && ~ischar(arg6)
    PosScalar = 'N';
    if isscalar(arg2)
        if PositiveScalars(arg2) % Check is scalar is positive. Otherwise, we have an error
           PosScalar = '2'; 
        end
    elseif isscalar(arg3)
        if PositiveScalars(arg3) % Check is scalar is positive. Otherwise, we have an error
           PosScalar = '3';
        end
    elseif isscalar(arg4)
        if PositiveScalars(arg4) % Check is scalar is positive. Otherwise, we have an error
           PosScalar = '4';
        end
    elseif isscalar(arg5)
        if PositiveScalars(arg5) % Check is scalar is positive. Otherwise, we have an error
           PosScalar = '5';
        end
    end
    
    switch PosScalar
        case '2'
            if CheckNextScalar(arg3) % Check whether next arguement is a postive scalar. Otherwise, we have an error
                rhoMu = arg2;
                rhoSigma = arg3;
                [SupportSet,MeanCenter,VarianceCenter] = CheckRemaining3Args(arg4,arg5,arg6);                
            end
        case '3'
            if CheckNextScalar(arg4) % Check whether next arguement is a postive scalar. Otherwise, we have an error
                rhoMu = arg3;
                rhoSigma = arg4;
                [SupportSet,MeanCenter,VarianceCenter] = CheckRemaining3Args(arg2,arg5,arg6);                
            end
        case '4'
            if CheckNextScalar(arg5) % Check whether next arguement is a postive scalar. Otherwise, we have an error
                rhoMu = arg4;
                rhoSigma = arg5;
                [SupportSet,MeanCenter,VarianceCenter] = CheckRemaining3Args(arg2,arg3,arg6);
            end
        case '5'
            if CheckNextScalar(arg6) % Check whether next arguement is a postive scalar. Otherwise, we have an error
                rhoMu = arg5;
                rhoSigma = arg6;
                [SupportSet,MeanCenter,VarianceCenter] = CheckRemaining3Args(arg2,arg3,arg4);
            end
        otherwise
            error('We must have exactly two scalars as input')
    end
else
    error('At least one of the arguments must be a real number')
end

% Error if the dimension of Sigma is not equal to the rows of the
% SupportSet
if length(VarianceCenter) ~= size(SupportSet,1)
    error('The dimension of the input matrix must be equal to the number of rows of the cardinality set');
end

% Error if the dimension of ObjFunc is not equal to the columns of the
% SupportSet
if length(ObjFunc) ~= size(SupportSet,2)
    error('The dimension of the cost function must be equal to size of the support set');
end

% Converting MeanCenter into column vector if necessary
if ~isempty(MeanCenter) && size(MeanCenter,1) == 1
    MeanCenter = MeanCenter';
end
% 
% % Checking whether MeanCenter is a probability measure
% if sum(MeanCenter < 0) ~= 0 || abs(sum(MeanCenter) - 1) > tol 
%     error('The mean center must be a probability vector')    
% end

end

function out = CheckNextScalar(arg1)
tol = 1e-5;

if isscalar(arg1) && arg1 > tol
    out = true;
else
    error('Either next input is not scalar or it is a negative scalar')
end

end

function out = PositiveScalars(arg1)

tol = 1e-5;

if arg1 < -tol
    error('All scalars must be positive')
else
    out = true;
end
end

function [SupportSet,MeanCenter,VarianceCenter] = CheckRemaining3Args(arg1,arg2,arg3)

SupportSet = arg3;
[MeanCenter,VarianceCenter] =  CheckRemaining2Args(arg1,arg2);

if ismatrix(SupportSet) && ~ischar(SupportSet)
    if size(SupportSet,1) ~= length(MeanCenter)
        error('The dimension of each element in the support set must be the same as that of the mean center')
    end
else
    error('The support set must be a row or column vector of numbers')
end

end

function [MeanCenter,VarianceCenter] = CheckRemaining2Args(arg1,arg2)

    tol = 1e-5;
    if isvector(arg1) && ~ischar(arg1) || isvector(arg2) && ~ischar(arg2)
        if isvector(arg1)
            MeanCenter = arg1;
            if CheckPSDandSymmetric(arg2)
                VarianceCenter = arg2;
            end
        else
            MeanCenter = arg2;
            if CheckPSDandSymmetric(arg1)
                VarianceCenter = arg1;
            end
        end
    else
        error('At least one argument must be a row or column vector of numbers')
    end

end

function out = CheckPSDandSymmetric(arg1)

tol = 1e-5;

if norm(arg1 - arg1') < tol && min(real(eig(arg1))) > -tol && (ismatrix(arg1) && (size(arg1,1) == size(arg1,2)))
    out = true;
else
    error('The matrix must be square, symmetric and PSD')
end

end




