function value = RemainingIterations(NumberOfNestedLoops,indicesAndLength,TailIterations,barHandle)

if find(indicesAndLength(:,1) < 0)
    close(barHandle);
    error('All indices must be larger than zero');
end

if NumberOfNestedLoops ~= length(indicesAndLength)
    close(barHandle);
    error('The first argument must be equal to the dimension of the second argument');
end

tempValue = 0;

for j = 1:NumberOfNestedLoops  
    if j~= NumberOfNestedLoops
        tempValue = tempValue + prod(indicesAndLength(j+1:end,2))*(indicesAndLength(j,1)-1)*TailIterations;
    else
        tempValue = tempValue + indicesAndLength(end,1)*TailIterations;
    end
end

value = double(tempValue);

end