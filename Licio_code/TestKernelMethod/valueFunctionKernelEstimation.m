function out = valueFunctionKernelEstimation(horizon,Data,statePart,safeSet,sigma,lambda,eta)

m = size(Data,1);

PsiEval = zeros(m,m);

ValueFunc = zeros(m,horizon+1);
ValueFunc(:,end) = initializeValueFunc(statePart,safeSet);

for i = 1:m
    currentState = Data(i,2);
    PsiEval(:,i) = evalConditionalMeanEmbedding(currentState,statePart,Data,sigma);
end

A = PsiEval + lambda*m*eye(m);
b = PsiEval;

temp = A\b;

for i=horizon:-1:1
    tempVal = (ValueFunc(:,i+1)'*temp)';
    ValueFunc(:,i) = eta*min(1,max(0,tempVal));
end

out = ValueFunc(:,1);


end