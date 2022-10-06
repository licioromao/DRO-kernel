function [] = PrintOuterLoopIterations(TotalIterations,CurrentIteration,TimeCurrentIteration)

fprintf('\n')

fprintf('Total of iterations: %d \n',TotalIterations);

fprintf('\n')
disp('|===================================================================================================|')
NString = 100;

perc_iterates = CurrentIteration/TotalIterations;

progress_string = '';
for i =1:floor(perc_iterates*NString)
    progress_string = [progress_string,'#'];
end

disp(progress_string);
disp(['|========================================= ',num2str(floor(perc_iterates*100)),'% completed ===========================================|']);

%display progress per cent
steps_remaining = TotalIterations - CurrentIteration;
minutes = floor(TimeCurrentIteration * steps_remaining / 60);
seconds = rem(floor(TimeCurrentIteration *  steps_remaining), 60);

if CurrentIteration > 1
    disp(' ');
    
    disp(['  Estimated remaining time: ', num2str(minutes), '(min) ', num2str(seconds), '(sec) ']);
    
    fprintf('\n')
    fprintf('\n')
end

end

