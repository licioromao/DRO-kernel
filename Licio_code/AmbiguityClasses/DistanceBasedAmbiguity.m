classdef DistanceBasedAmbiguity < Ambiguity
    
    % This class is an abstract subclass of the abstract class ambiguity and
    % will implement the abstract method SolveOptimization using a distance in
    % the space of probability measures. To this end, in addition to the
    % parameters of the class ambiguity, this abstract class possesses two
    % other parameters (center_ball and radius_ball) that defines the center 
    % and radius of the ambiguity set
   
    properties (Access = protected)
        center_ball % Probability vector that is the center of the ambiguity set 
        radius_ball % Radius of the ball
    end
    
    methods
        function obj = DistanceBasedAmbiguity(objective_cost,radius,center)
                    
            [~,~,~] = verify_arg_distance_based_class(objective_cost,radius,center);

            obj = obj@Ambiguity(objective_cost); % Call the constructor of the parent class
            obj.center_ball = center; % Initialize the center of the ambiguity set
            obj.radius_ball = radius; % radius of the ambiguity set
        end
        
        function out = get_values(obj)
            % Returns the private properties of this object
            out.objective_cost = obj.objective_cost; % This is a property of the parent class
            out.center_ball = obj.center_ball;
            out.radius_ball = obj.radius_ball;
        end
        
        function obj = set_values(obj,arg1,arg2,arg3)
            
            % Sets the value of the private properties the class. Note that
            % objective_cost is inhereted from the parent class
            [radius,center,objective_cost] = verify_arg_distance_based_class(arg1,arg2,arg3);
            
            obj.objective_cost = objective_cost;
            obj.center_ball = center; % Initialize the center of the ambiguity set
            obj.radius_ball = radius; % radius of the ambiguity set
        end
    end       
end

function [radius,center,objective_cost] = verify_arg_distance_based_class(arg1,arg2,arg3)

% As this abstract subclass contains three parameters (one
% inherited from its parent class) and two others that define
% the ambiguity set, we need to override the function
% verifyArg.

tol = 1e-5;

% Checking whether there exists a scalar as input parameter
if (isscalar(arg1)  && ~ischar(arg1))||(isscalar(arg2)  && ~ischar(arg2))||(isscalar(arg3)  && ~ischar(arg3))
    
    if isscalar(arg1)
        if arg1 < -tol
            error('The scalar parameter must be non-negative'); % an error if the scalar is negative (up to some tolerance)
        else
            radius = arg1; % set the radius as the positive scalar parameterr
        end
        SCALAR = '1';
    elseif isscalar(arg2)
        if arg2 < -tol
            error('The scalar parameter must be non-negative'); % an error if the scalar is negative (up to some tolerance)
        else
            radius = arg2; % set the radius as the positive scalar parameterr
        end
        SCALAR = '2';
    else
        if arg3 < -tol 
            error('The scalar parameter must be non-negative'); % an error if the scalar is negative (up to some tolerance)
        else
            radius = arg3; % set the radius as the positive scalar parameterr
        end
        SCALAR = '3';
    end
else
    error('One of the arguments must be a positive scalar number') % an error the radius is not passed as input
end

% Checking whether there exists a probability vector as input parameter
switch SCALAR
    case '1' % if scalar in the first argument, check arg2 and arg3
        [center,objective_cost] = verify_two_arg_distance_based_class(arg2,arg3); % Checking and saving second and thrid arguments
    case '2'
        [center,objective_cost] = verify_two_arg_distance_based_class(arg1,arg3);  % Checking and saving the first and third arguments
    case '3'
        [center,objective_cost] = verify_two_arg_distance_based_class(arg1,arg2);  % Checking and saving the first and second arguments
    otherwise
        error('More arguments than necessary to initialize this object');
end

if length(center) ~= length(objective_cost) && ~isempty(objective_cost) && ~isempty(center)
    error('Private variables c and q must be vectors of the same length');
end

end

function [center,objective_cost] = verify_two_arg_distance_based_class(arg1,arg2)
% This function tests whether one of the arguments is a
% probability measure. This is used to set values for the
% private parameters of the object

tol = 1e-5;

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
if (isempty(find(arg1 < 0)) && (abs(sum(arg1) - 1) <= tol)) || (isempty(find(arg2 < 0)) && (abs(sum(arg2) - 1) <= tol))
    if isempty(find(arg1 < 0)) && sum(arg1) == 1
        center = arg1;
        objective_cost = arg2;
    else
        center = arg2;
        objective_cost = arg1;
    end
else
    error('One of the vectors must be a probability measure');
end

end