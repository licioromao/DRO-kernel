function out = generatePartition(No_partition,SafeSet)

% Generates a lattice partition of the state space, taking as limit the
% safe set
%
% Input: No_Partition - a vector containing the number of points at each
%                       state of the model
%
%        SafeSet - A 3x2 matrix with the limits of the safe set in each
%                  dimension
%
% Output: out - is a structure with four fields as defined below:
%               out.X1, out.X2, out.X3 -- multi-dimensional matrices of size No_Partition(1)
%                                         x No_Partition(2) x No_Partition(3) containing the points
%                                         of the generated lattice
%               out.grid_x -- elements of the grid organized in a (No_Partition(1)
%                                         x No_Partition(2) x
%                                         No_Partition(3)) x 3 matrix
% The field out.grid_x is used to assign labels to the discrete states of
% the model.

grid_x1 = linspace(SafeSet(1,1),SafeSet(1,2),No_partition(1));
grid_x2 = linspace(SafeSet(2,1),SafeSet(2,2),No_partition(2));
grid_x3 = linspace(SafeSet(3,1),SafeSet(3,2),No_partition(3));

[X1,X2,X3] = meshgrid(grid_x1,grid_x2,grid_x3);

Nx1 = size(X1,2);
Nx2 = size(X1,1);
Nx3 = size(X1,3);

grid_x = zeros(Nx1*Nx2*Nx3,3);

for i1 = 1:Nx1
    for i2 = 1:Nx2
        for i3 = 1:Nx3
            x = [X1(i2,i1,i3);X2(i2,i1,i3);X3(i2,i1,i3)];
            index = (i1-1)*Nx2*Nx3 + (i2-1)*Nx3 + i3;
            grid_x(index,:) = x;
        end
    end
end

out.X1 = X1;
out.X2 = X2;
out.X3 = X3;
out.grid_x = grid_x;

end