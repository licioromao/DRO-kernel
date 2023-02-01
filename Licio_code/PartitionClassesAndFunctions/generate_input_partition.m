function input_partition = generate_input_partition(input,type_vector_field)
% This function generates a grid of size Nu^2 x 2 that contains any
% possible combination of discrete control inputs. The number of column is
% equal to the number of patches, which in this example is equal to 2.

switch type_vector_field
    case 'Fishery'
        number_of_inputs = length(input);
        input_partition = zeros(number_of_inputs^2,2); % Initializing the output
        
        for i=1:number_of_inputs
            for j =1:number_of_inputs
                input_partition(j+number_of_inputs*(i-1),1) = input(i);
                input_partition(j+number_of_inputs*(i-1),2) = input(j);
            end
        end
        
    case 'TCL'
        input_partition = [0;1];
    case 'ChainInt'
        input_partition = linspace(0,1,10)';
    case 'CarPole'
        input_partition = linspace(-10,10,5)';
    case 'CarPoleNL'
        input_partition = linspace(-10,10,5)';
    otherwise
        NoImplement();

end