function [] = PrintInnerLoop(TotalIteration,CurrentIteration,TimeCurrentIteration,CurrentAmbiguity,OuterLoopInfo)

clc;
StringAmbiguity = OuterLoopInfo.StringAmbiguity;
TotalNumberAmbiguity = length(StringAmbiguity);

fprintf('Executing vector field %s with the following ambiguity sets:\n',OuterLoopInfo.TypeVectorField)

for i=1:length(StringAmbiguity)
    fprintf('%d. %s\n',i,StringAmbiguity{i}.Name)
end

PrintOuterLoopIterations(OuterLoopInfo.TotalIteration,OuterLoopInfo.CurrentIteration,OuterLoopInfo.TimeIteration);

key = false;
key2 = false;
i = 1;

while ~key && i <= TotalNumberAmbiguity
    if ~strcmp(CurrentAmbiguity,StringAmbiguity{i}.Name)
        if ~key2
            disp('(Inner loop) Ambiguity sets completed:')
            key2 = true;
        end
        fprintf('%d. %s\n',i -1,StringAmbiguity{i}.Name)
    else
        key = true;
    end
    i = i + 1;
end

fprintf('\n')
fprintf('Computing value function of %s \n', CurrentAmbiguity);
fprintf('(Inner loop) Total of iterations: %d \n',TotalIteration);

fprintf('\n')
disp('|================================================|')
NString = 50;

perc_iterates = CurrentIteration/TotalIteration;

progress_string = '';
for i =1:floor(perc_iterates*NString)
    progress_string = [progress_string,'#'];
end

disp(progress_string);
disp(['|================ ',num2str(floor(perc_iterates*100)),'% completed =================|']);

%display progress per cent
steps_remaining = TotalIteration - CurrentIteration;
minutes = floor(TimeCurrentIteration * steps_remaining / 60);
seconds = rem(floor(TimeCurrentIteration *  steps_remaining), 60);
disp(' ');

disp(['  Estimated remaining time: ', num2str(minutes), '(min) ', num2str(seconds), '(sec) ']);

fprintf('\n')
fprintf('\n')


end

