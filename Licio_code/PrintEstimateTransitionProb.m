function [] = PrintEstimateTransitionProb(TotalIterations,CurrentIteration,TimeCurrentIteration)

clc;

disp('Generating the transition probability.');

fprintf('Number of states-action pair: %d \n', TotalIterations);

PrintOuterLoopIterations(TotalIterations,CurrentIteration,TimeCurrentIteration);

end

