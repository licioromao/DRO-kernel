function final_value_func = initialise_value_func(state_partition_grid,Q)

number_of_points = size(state_partition_grid,1);
final_value_func = zeros(number_of_points,1);

for i =1:number_of_points
    final_value_func(i) = -state_partition_grid(i,:)*Q*state_partition_grid(i,:)';
end

end