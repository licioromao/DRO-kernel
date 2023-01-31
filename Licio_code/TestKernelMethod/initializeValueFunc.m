function ValueFunc = initializeValueFunc(StatePart,safeSet)


grid_X = StatePart.getValues.Partition.grid_x;
m = length(grid_X);

ValueFunc = zeros(m,1);

for i=1:m
    if safeSet(1) < grid_X(i) && grid_X(i) < safeSet(2)
        ValueFunc(i) = 1;
    end
end

end