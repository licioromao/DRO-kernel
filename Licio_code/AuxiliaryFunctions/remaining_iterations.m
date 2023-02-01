function value = remaining_iterations(number_of_nested_loops,indices_and_length,...
                                            tail_iterations,bar_handle)

if find(indices_and_length(:,1) < 0)
    close(bar_handle);
    error('All indices must be larger than zero');
end

if number_of_nested_loops ~= size(indices_and_length,1)
    close(bar_handle);
    error('The first argument must be equal to the dimension of the second argument');
end

temp_value = 0;

for j = 1:number_of_nested_loops  
    if j~= number_of_nested_loops
        temp_value = temp_value + prod(indices_and_length(j+1:end,2))*(indices_and_length(j,1)-1)*tail_iterations;
    else
        temp_value = temp_value + indices_and_length(end,1)*tail_iterations;
    end
end

value = double(temp_value);

end