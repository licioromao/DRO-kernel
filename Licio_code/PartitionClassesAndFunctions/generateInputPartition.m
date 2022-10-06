function out = generateInputPartition(Input,TypeOfVectorField)
% This function generates a grid of size Nu^2 x 2 that contains any
% possible combination of discrete control inputs. The number of column is
% equal to the number of patches, which in this example is equal to 2.

switch TypeOfVectorField
    case 'Fishery'
        Nu = length(Input);
        out = zeros(Nu^2,2); % Initializing the output
        
        for i=1:Nu
            for j =1:Nu
                out(j+Nu*(i-1),1) = Input(i);
                out(j+Nu*(i-1),2) = Input(j);
            end
        end
        
    case 'TCL'
        out = [0;1];
    case 'ChainInt'
        out = linspace(0,1,10)';
    case 'CarPole'
        out = linspace(-10,10,5)';
    case 'CarPoleNL'
        out = linspace(-10,10,5)';
    otherwise
        NoImplement();

end