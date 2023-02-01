function [] = print_estimate_transition_prob(total_iterations,current_iteration...
                                                ,time_current_iteration)

clc;

disp('Generating the transition probability.');

fprintf('Number of states-action pair: %d \n', total_iterations);

print_outer_loop_iterations(total_iterations,current_iteration,time_current_iteration);

end

