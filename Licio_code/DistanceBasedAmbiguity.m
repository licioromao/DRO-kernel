classdef DistanceBasedAmbiguity < Ambiguity
    
    % This class is an abstract subclass of the abstract class ambiguity and
    % will implement the abstract method SolveOptimization using a distance in
    % the space of probability measures. To this end, in addition to the
    % parameters of the class ambiguity, this abstract class possess two
    % other parameters (p and epsilon) that defines the center and radius
    % of the ambiguity set
   
    properties (Access = protected)
        q % Probability vector that is the center of the ambiguity set 
        epsilon % Radius of the ball
    end
    
    methods
        function obj = DistanceBasedAmbiguity(arg1,arg2,arg3)
                    
            [Epsilon,ProbMeas,ObjFunc] = verifyArgDistanceBased(arg1,arg2,arg3);
            obj = obj@Ambiguity(ObjFunc); %Call the constructor of the parent class
            
            obj.q = ProbMeas; % Initialize the center of the ambiguity set
            obj.epsilon = Epsilon; % radius of the ambiguity set
        end
        
        function out = getValues(obj)
            % Returns the private properties of this object
            out.c = obj.c; % This is a property of the parent class
            out.q = obj.q;
            out.epsilon = obj.epsilon;
        end
        
        function obj = setValues(obj,arg1,arg2,arg3)
            
            % Sets the value of the private properties the class. Note that
            % c is inhereted from the parent class
            
            [Epsilon,ProbMeas,ObjFunc] = verifyArgDistanceBased(arg1,arg2,arg3);
            
            obj.c = ObjFunc;
            obj.q = ProbMeas; % Initialize the center of the ambiguity set
            obj.epsilon = Epsilon; % radius of the ambiguity set
        end
    end       
end

function [epsilon,ProbMeas,ObjFunc] = verifyArgDistanceBased(arg1,arg2,arg3)

% As this abstract subclass contains three parameters (one
% inherited from its parent class) and two others that define
% the ambiguity set, we need to override the function
% verifyArg.

tol = 1e-5;

% Checking whether there exists a scalar as input parameter
if (isscalar(arg1)  && ~ischar(arg1))||(isscalar(arg2)  && ~ischar(arg2))||(isscalar(arg3)  && ~ischar(arg3))
    scalar = 'N';
    if isscalar(arg1)
        if arg1 < -tol
            error('The scalar parameter must be non-negative'); % an error if the scalar is negative (up to some tolerance)
        else
            epsilon = arg1; % set the radius as the positive scalar parameterr
        end
        scalar = '1';
    elseif isscalar(arg2)
        if arg2 < -tol
            error('The scalar parameter must be non-negative'); % an error if the scalar is negative (up to some tolerance)
        else
            epsilon = arg2; % set the radius as the positive scalar parameterr
        end
        scalar = '2';
    else
        if arg3 < -tol 
            error('The scalar parameter must be non-negative'); % an error if the scalar is negative (up to some tolerance)
        else
            epsilon = arg3; % set the radius as the positive scalar parameterr
        end
        scalar = '3';
    end
else
    error('One of the arguments must be a positive scalar number') % an error the radius is not passed as input
end

% Checking whether there exists a probability vector as input parameter
switch scalar
    case '1' % if scalar is the first argument, check arg2 and arg3
        [ProbMeas,ObjFunc] = VerifyTwoArg(arg2,arg3); % Checking and saving second and thrid arguments
    case '2'
        [ProbMeas,ObjFunc] = VerifyTwoArg(arg1,arg3);  % Checking and saving the first and third arguments
    case '3'
        [ProbMeas,ObjFunc] = VerifyTwoArg(arg1,arg2);  % Checking and saving the first and second arguments
    otherwise
        error('More arguments than necessary to initialize this object');
end

if length(ProbMeas) ~= length(ObjFunc) && ~isempty(ObjFunc) && ~isempty(ProbMeas)
    error('Private variables c and q must be vectors of the same length');
end
end

function [ProbMeas,ObjFunc] = VerifyTwoArg(arg1,arg2)
            % This function tests whether one of the arguments is a
            % probability measure. This is used to set values for the
            % private parameters of the object
            
            if ~(isvector(arg1) && isvector(arg2) && ~ischar(arg1) && ~ischar(arg2))
                error('Two of the inputs must be column or row vectors of real numbers')
            end
                       
            % Transform the inputs into column vectors
            if size(arg1,1) == 1
                arg1 = arg1';
            end
            
            if size(arg2,1) == 1
                arg2 = arg2';
            end
            
            % Checking whether there exists a prob measure amongst the
            % inputs; otherwise, ouput an error
            if (isempty(find(arg1 < 0)) && sum(arg1) == 1) || (isempty(find(arg2 < 0)) && sum(arg2) == 1) 
                if isempty(find(arg1 < 0)) && sum(arg1) == 1
                    ProbMeas = arg1;
                    ObjFunc = arg2;
                else
                    ProbMeas = arg2;
                    ObjFunc = arg1;
                end
            else
                error('One of the vectors must be a probability measure');
            end
        end