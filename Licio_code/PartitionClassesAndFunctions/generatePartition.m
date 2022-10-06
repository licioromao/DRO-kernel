function out = generatePartition(No_partition,SafeSet,typeOfVectorField)

% Generates a lattice partition of the state space, taking as limit the
% safe set
%
% Input: No_Partition - a vector containing the number of points at each
%                       state of the model
%
%        SafeSet - A 3x2 matrix with the limits of the safe set in each
%                  dimension
%        typeOfVectorField - This is a string containing the type of vector field that are implemented. Currently, we support only Fishery and TCL vector fields.
%
% Output: out - is a structure with different fields depending on the type of vector field parameter. If this is equal to 'Fishesry' then we have that:
%               out.X1, out.X2, out.X3 -- multi-dimensional matrices of size No_Partition(1)
%                                         x No_Partition(2) x No_Partition(3) containing the points
%                                         of the generated lattice
%               out.grid_x -- elements of the grid organized in a (No_Partition(1)
%                                         x No_Partition(2) x
%                                         No_Partition(3)) x 3 matrix
%                If the type of vector field parameter is equal to 'TCL', then the output will be given by:
%                out.X -- Nx1 points containing the grid of the one-dimensional space
%
%                out.grid_x -- This coincides with out.X in this case, since the space has one dimension.
%
%
%
%
% The field out.grid_x is used to assign labels to the discrete states of
% the model.
switch typeOfVectorField
    case '1D'
        ub = max(SafeSet) + 1;
        lb = min(SafeSet) - 1;
        h = (ub-lb)/(No_partition-1);
        
        tempOut = lb:h:ub;
        
        out.X = tempOut';
        out.grid_x = out.X;
        
    case '2D'
        
        grid_x1 = linspace(SafeSet(1,1),SafeSet(1,2),No_partition(1));
        grid_x2 = linspace(SafeSet(2,1),SafeSet(2,2),No_partition(2));
        
        [X1,X2] = meshgrid(grid_x1,grid_x2);
        
        Nx1 = size(X1,2);
        Nx2 = size(X1,1);
        
        grid_x = zeros(Nx1*Nx2,2);
        
        for i1 = 1:Nx1
            for i2 = 1:Nx2
                x = [X1(i2,i1);X2(i2,i1)];
                index = RemainingIterations(2,[[i1;i2],[Nx1;Nx2]],1,[]);
                grid_x(index,:) = x;
            end
        end
        
        out.X1 = X1;
        out.X2 = X2;
        out.grid_x = grid_x;
        
    case '3D'
        
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
     
    case '4D'
        
        grid_x1 = linspace(SafeSet(1,1),SafeSet(1,2),No_partition(1));
        grid_x2 = linspace(SafeSet(2,1),SafeSet(2,2),No_partition(2));
        grid_x3 = linspace(SafeSet(3,1),SafeSet(3,2),No_partition(3));
        grid_x4 = linspace(SafeSet(4,1),SafeSet(4,2),No_partition(4));
        
        [X1,X2,X3,X4] = ndgrid(grid_x1,grid_x2,grid_x3,grid_x4);
        
        
        Nx1 = size(X1,1);
        Nx2 = size(X1,2);
        Nx3 = size(X1,3);
        Nx4 = size(X1,4);
        
        grid_x = zeros(Nx1*Nx2*Nx3*Nx4,4);
        
        for i1 = 1:Nx1
            for i2 = 1:Nx2
                for i3 = 1:Nx3
                    for i4=1:Nx4
                        x = [X1(i1,i2,i3,i4);X2(i1,i2,i3,i4);X3(i1,i2,i3,i4);X4(i1,i2,i3,i4)];
                        index = RemainingIterations(4,[[i1;i2;i3;i4],[Nx1;Nx2;Nx3;Nx4]],1,[]);
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
        NotImplemented();
end
        
end
