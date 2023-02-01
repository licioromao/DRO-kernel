function [] = print_inner_loop(total_iteration,current_iteration,time_current_iteration...
                                ,current_ambiguity,outer_loop_info)

clc;
string_ambiguity = outer_loop_info.string_ambiguity;
total_number_ambiguity = length(string_ambiguity);

fprintf('Executing vector field %s with the following ambiguity sets:\n',outer_loop_info.type_vector_field)

for i=1:length(string_ambiguity)
    fprintf('%d. %s\n',i,string_ambiguity{i}.name)
end

print_outer_loop_iterations(outer_loop_info.total_iteration,...
                                outer_loop_info.current_iteration,...
                                    outer_loop_info.time_iteration);

key = false;
key2 = false;
i = 1;

while ~key && i <= total_number_ambiguity
    if ~strcmp(current_ambiguity,string_ambiguity{i}.name)
        if ~key2
            disp('(Inner loop) Ambiguity sets completed:')
            key2 = true;
        end
        fprintf('%d. %s\n',i -1,string_ambiguity{i}.name)
    else
        key = true;
    end
    i = i + 1;
end

fprintf('\n')
fprintf('Computing value function of %s \n', current_ambiguity);
fprintf('(Inner loop) Total of iterations: %d \n',total_iteration);

fprintf('\n')
disp('|================================================|')
NString = 50;

perc_iterates = current_iteration/total_iteration;

progress_string = '';
for i =1:floor(perc_iterates*NString)
    progress_string = [progress_string,'#'];
end

disp(progress_string);
disp(['|================ ',num2str(floor(perc_iterates*100)),'% completed =================|']);

%display progress per cent
steps_remaining = total_iteration - current_iteration;
minutes = floor(time_current_iteration * steps_remaining / 60);
seconds = rem(floor(time_current_iteration *  steps_remaining), 60);
disp(' ');

disp(['  Estimated remaining time: ', num2str(minutes), '(min) ', num2str(seconds), '(sec) ']);

fprintf('\n')
fprintf('\n')


end

