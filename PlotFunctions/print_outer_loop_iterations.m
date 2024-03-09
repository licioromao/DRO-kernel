function [] = print_outer_loop_iterations(total_iterations,current_iteration,...
                                            time_current_iteration)

fprintf('\n')

fprintf('Total of iterations: %d \n',total_iterations);

fprintf('\n')
disp(['|======================================================================='...
                        '============================|'])
NString = 100;

perc_iterates = current_iteration/total_iterations;

progress_string = '';
for i =1:floor(perc_iterates*NString)
    progress_string = [progress_string,'#'];
end

disp(progress_string);
disp(['|========================================= ',...
    num2str(floor(perc_iterates*100)),'% completed =====================' ...
                                                    '======================|']);

%display progress per cent
steps_remaining = total_iterations - current_iteration;
minutes = floor(time_current_iteration * steps_remaining / 60);
seconds = rem(floor(time_current_iteration *  steps_remaining), 60);

if current_iteration > 1
    disp(' ');
    
    disp(['  Estimated remaining time: ', num2str(minutes), '(min) ', ...
                                                    num2str(seconds), '(sec) ']);
    
    fprintf('\n')
    fprintf('\n')
end

end

