% This file defines a class that is used to compute transition
% probabilities for a MDP abstraction of continuous-time models.
%
% The idea for creating this class is motivated by the need to simulate
% transitions from a given state-action pair and produce a sampling-based
% estimate of the transition probabilities to the next state. This class
% contains public and private field, as well as several methods to
% manupulate each of these fields.
%
%
%   Public fields: List -- contains the set of allowable names for the list
%                          as a cell array of strings
%                  lenghList -- number of allowable names
%
%   Private fields: CurrentList -- Stores the elements that are present in
%                                  the list. All its elements must be an
%                                  element of the allowable public field
%                   Value -- This is a number associated with an element
%                           of the private list. In our application, this value
%                           corresponding to the number of transitions to
%                           each element of the MDP states.
%
%   Intuition to create this class: All possible state-action pairs of the
%   MDP will be present in the public field List and the private fields
%   will contain the probability measure of transitioning to the next state
%   of the MDP.

classdef HashTable
    
    % Public parameters
    properties
        List 
        lengthList % length of the list
    end
    
    % Private parameters
    properties (Access = protected)
        CurrentList
        Value % Contains the value associated with each element of the table
    end
    
    methods (Static, Access = private)
        
        % The use of this private function is to test the arguments that
        % are used to change the values of this hash table
        
        function [] = testSetValueArg(arg1,arg2)
            
            % This private function checks whether arg1 and arg2 have
            % compatible dimension. We require arg1 to be a cell-array with
            % entries given by strings, and arg2 is a column or row vector
            % whose size is equal to the length of arg1.
            
            out(1) = false; 
            out(2) = true;
            
            % Checking whether arg1 has the correct property
            if iscell(arg1)
                n = length(arg1); % length of arg1
                b = true; % temporary boolean variable
                for i = 1:n
                    b = b && ischar(arg1{i});
                end % If all elements of arg1 are strings then b == true in the end of this for loop
                
                % Update the first output if arg1 passes the above test
                if b
                    out(1) = true;
                end
            end
            
            if out(1) % If out(1) is true, then check whether arg2 has the same length
                if length(arg1) ~= length(arg2)
                    out(2) = false;
                end
            end
            
            % Output an error is one is found
            if ~out(1)
                error('The first argument is either not a cell, or it is not a cell whose entries are all strings');
            end
            
            if ~out(2)
                error('The length of cell must coincide with the length of the second argument');
            end
            
        end
    end
    
    methods (Access = private)
        function out = checkMembershipPrivateList(obj,arg1)
            
            % This function is similar to the other private function
            % checkMembershipMainList in that it verifies whether arg1 is
            % an element of the private list, returning an error if arg1 is
            % duplicated.
            
            a = sum(strcmp(obj.CurrentList,arg1));
            if a == 0
                out = false;
            elseif a == 1
                out = true;
            else
                error('There must be an error with the data structure. We cannot store an element twice.')
            end
        end
        
        function out = checkMembershipMainList(obj,arg1)
            
            % This function checks if the string arg1 is an element in the
            % main list and returns an error if the element is present
            % twice.
            
            a = sum(strcmp(obj.List,arg1)); 
            if a == 0
                out = false;
            elseif a == 1
                out = true;
            else
                error('There must be an error with the data structure. The field List cannot contain the same element twice.') % Outputs an error if an element is present in the list twice
            end
        end
        
        function [] = testStructure(obj)
            
            % This private function checks consistency of the list. The
            % current version verifies whether the private and public fields have
            % compatible lengths; it also checks whether all elements in
            % the field CurrentList is present in the field List.
            
            % Compatible private and public fields
            if length(obj.CurrentList) ~= length(obj.Value) || length(obj.List) ~= obj.lengthList
                error('There is an error with this Hastable. Please create another object and delete this one.')
            end
            
            % Whenever the public list is non-empty, the test below makes
            % sure all elements in the private list is an element of the
            % public one.
            if ~isempty(obj.List)
                L = length(obj.CurrentList);
                for i =1:L
                    a = obj.checkMembershipMainList(obj.CurrentList{i});
                    if ~a
                        error('The internal list contains an element that does not belong to the main List. Please create another object and delete this one');
                    end
                end
            end
            
        end
        
    end
    
    
    methods
        function obj = HashTable(arg1)
            % Constructor function of the class
            if isempty(arg1)
                obj.List = {};
                obj.CurrentList = {};
                obj.Value = [];
                obj.lengthList = 0;
            else
                obj.CurrentList = arg1;
                obj.Value = zeros(length(obj.CurrentList),1) ;
                obj.List = {};
                obj.lengthList = 0;
            end
            
            obj.testStructure();
        end
        
        function obj = setValues(obj,arg1,arg2)
            
            % This function attempts to append a new element arg1 (if not yet
            % present) with value equal to arg2 to the internal list. In doing so, it checks
            % whether these have compatible dimensions and whether this is
            % compatible with the values in the field List.
            
             obj.testSetValueArg(arg1,arg2); % Checking whether arg1 and arg2 are compatible
             
            
            % This checks whether arg1 is a roww or column vector cell
            % array. If is a row vector, we transpose it. If not a row or
            % collumn vector array, we output an error.
            if min(size(arg1)) ~= 1
                error('The input must be a row or column cell array') % Error in case arg1 is not a row and column array
            elseif size(arg1,1) == 1 % If a row array, then we transpose it
                tempArg1 = arg1';
            else
                tempArg1 = arg1;
            end
            
            if ~isempty(obj.List)
                for i=1:length(tempArg1)
                    a = obj.checkMembershipMainList(tempArg1{i});
                    if ~a
                        error('You cannot set the internal list to these strings, as there is an element of it that does not belong to the main list'); % Unable to append arg1 if not present in the main list
                    end
                end
            else
                error('We must initialize the public field List before modifying the internal list') % Error if list of allowable names are not defined
            end
            
            obj.CurrentList = arg1;
            obj.Value = arg2;
            
        end
        
        function out = getValues(obj)
            
            % Returns the values of the private variables for manipulations
            
            out.CurrentList = obj.CurrentList;
            out.Value = obj.Value;
        end
        
        function obj = appendString(obj,arg1)
            
            % Appends an element to the private list. 
            
            % Verifies whether arg1 is a string
            if ~ischar(arg1)
                error('The argument of this function must be a string');
            end
            
            % Checks whether arg1 is an allowable string
            if ~isempty(obj.List)
                a = obj.checkMembershipMainList(arg1);
                if ~a
                    error('You cannot add elements that are not in the main List.')
                end
            end
            
            % Checks whether arg1 is already an element of the private list
            a = obj.checkMembershipPrivateList(arg1);
            
            if ~a % if not yet a member, create a new entry; otherwise, do nothing
                obj.CurrentList = cat(1,obj.CurrentList,arg1);
                obj.Value = [obj.Value;0]; % Initilize the corresponding value to zero
            end
            
            obj.testStructure(); % Test whether the modified list is admissible
        end
        
        function obj = removeString(obj,arg1)
            
            % Removes arg1 from the list if such an element exists
            
            a = obj.checkMembershipPrivateList(arg1); % Checks whether arg1 is an element of the list
            
            if a % If arg1 is an element, remove from the list
                index = find(strcmp(obj.CurrentList,arg1));
                obj.CurrentList(index) = [];
                obj.Value(index) = [];
            else % otherwise, do nothing and output an error
                warning('The argument is not an element of the list. No action has been taken')
            end
            
            obj.testStructure(); % Checks whether the resulting list is allowable
            
        end
        
        function obj = addValue(obj,arg1,arg2)
            
            % Adds the value arg2 to the private value associated with
            % arg1. First, we need to check whether arg1 is an element of
            % the private list.
            
            
           % Checking whether arg1 is a string
            if ~ischar(arg1)
                error('The first argument of this function must be a string');
            end
            
            a = obj.checkMembershipPrivateList(arg1); 
            
            if a % if arg1 is an element of the private list, add arg2 to its corresponding value
                index = find(strcmp(obj.CurrentList,arg1));
                obj.Value(index) = obj.Value(index) + arg2;
            else % otherwise, output an error
                error('Cannot add this value since %s is not an element of the list',arg1)
            end
            
            obj.testStructure(); % Checks whether the resulting list is admissible
        end
         
        function out = createProbMeasure(obj)
            
            % Returns the empirical transition probability to each element of the
            % list
            
            if isempty(obj.List)
                error('Please initialize the field List before calling this function'); % error if the public parameter List is empty
            else
                obj.testStructure(); % Checks whether we have a valid list            
                L = length(obj.Value);
                
                out = zeros(obj.lengthList,1); % Total number of possible transitions
                
                for i = 1:L
                    temp = obj.CurrentList{i};
                    
                    %Checks whether the list is consistent
                    obj.checkMembershipPrivateList(temp);
                    a = obj.checkMembershipMainList(temp); 
                    
                    if ~a
                        error('There is an error with this data structure. Do not trust the result of this function');
                    else % If not error is found, compute the corresponding probability measure
                        index = find(strcmp(obj.List,temp));
                        out(index) = obj.Value(i)/sum(obj.Value);
                    end
                end
                
                % Check whether we have a valid probability measure vector
                if sum(find(out < -1e-2)) > 0 || abs(sum(out) - 1) > 1e-2
                    error('The must be some error on the data structure. The resulting vector is not a probability distribution');
                end 
            end
            
        end
        
    end    
    
end

