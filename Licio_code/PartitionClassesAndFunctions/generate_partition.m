function out = generate_partition(grid,number_of_points)

% Generates a lattice partition of the state space, taking as limit the
% safe set
%
% Input: grid - A matrix with upper and lower limits for the grid
%        
%        number_of_points - a vector containing the number of points at each
%                       state of the model
%
%        safe_set - A matrix with the limits of the safe set in each dimension
%                  
%
% Output: out - is a structure with different fields depending on the dimension.
%               For instance, if dimension equal to three we have:
%
%               out.X1, out.X2, out.X3 -- multi-dimensional matrices of size number_of_points(1)
%                                         x number_of_points(2) x number_of_points(3) containing the points
%                                         of the generated lattice
%               out.grid_x -- elements of the grid organized in a (number_of_points(1)
%                                         x number_of_points(2) x
%                                         number_of_points(3)) x 3 matrix
%
% The field out.grid_x is used to assign labels to the discrete states of
% the model.

dimension = size(grid,1);

switch dimension
    case 1
        lb = grid(:,1);
        ub = grid(:,2);
        h = (ub-lb)/(number_of_points-1);
        
        temp_out = lb:h:ub;
        
        out.X = temp_out';
        out.grid_x = out.X;
        
    case 2
        
        grid_axis_1 = linspace(grid(1,1),grid(1,2),number_of_points(1));
        grid_axis_2 = linspace(grid(2,1),grid(2,2),number_of_points(2));
        
        [X1,X2] = meshgrid(grid_axis_1,grid_axis_2);
        
        grid_x = zeros(prod(number_of_points),2);
        
        for i1 = 1:number_of_points(1)
            for i2 = 1:number_of_points(2)
                x = [X1(i2,i1);X2(i2,i1)];
                index = remaining_iterations(2,[[i1;i2],number_of_points]...
                               ,1,[]); % Finds the number of remaining iterations
                grid_x(index,:) = x;
            end
        end
        
        out.X1 = X1;
        out.X2 = X2;
        out.grid_x = grid_x;
        
    case 3
        
        grid_axis_1 = linspace(grid(1,1),grid(1,2),number_of_points(1));
        grid_axis_2 = linspace(grid(2,1),grid(2,2),number_of_points(2));
        grid_axis_3 = linspace(grid(3,1),grid(3,2),number_of_points(3));
        [X1,X2,X3] = meshgrid(grid_axis_1,grid_axis_2,grid_axis_3);
        
        grid_x = zeros(prod(number_of_points),3);
        
        for i1 = 1:number_of_points(1)
            for i2 = 1:number_of_points(2)
                for i3 = 1:number_of_points(3)
                    x = [X1(i2,i1,i3);X2(i2,i1,i3);X3(i2,i1,i3)];
                    index = remaining_iterations(3,[[i1;i2;i3],number_of_points]...
                               ,1,[]); % Finds the number of remaining iterations
                    grid_x(index,:) = x;
                end
            end
        end
        
        out.X1 = X1;
        out.X2 = X2;
        out.X3 = X3;
        out.grid_x = grid_x;
     
    case 4
        
        grid_axis_1 = linspace(grid(1,1),grid(1,2),number_of_points(1));
        grid_axis_2 = linspace(grid(2,1),grid(2,2),number_of_points(2));
        grid_axis_3 = linspace(grid(3,1),grid(3,2),number_of_points(3));
        grid_axis_4 = linspace(grid(4,1),grid(4,2),number_of_points(4));
        
        [X1,X2,X3,X4] = ndgrid(grid_axis_1,grid_axis_2,grid_axis_3,grid_axis_4);
        
        grid_x = zeros(prod(number_of_points),4);
        
        for i1 = 1:number_of_points(1)
            for i2 = 1:number_of_points(2)
                for i3 = 1:number_of_points(3)
                    for i4=1:number_of_points(4)
                        x = [X1(i1,i2,i3,i4);X2(i1,i2,i3,i4);X3(i1,i2,i3,i4);X4(i1,i2,i3,i4)];
                        index = remaining_iterations(4,[[i1;i2;i3;i4],...
                                                        number_of_points],1,[]);

                        grid_x(index,:) = x;
                    end
                end
            end
        end
        
        out.X1 = X1;
        out.X2 = X2;
        out.X3 = X3;
        out.X4 = X4;
        out.grid_x = grid_x;
        
    otherwise
        not_implemented();
end
        
end
