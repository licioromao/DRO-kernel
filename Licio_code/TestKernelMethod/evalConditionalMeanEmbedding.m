function out = evalConditionalMeanEmbedding(currentState,StatePart,Data,sigma)

gamma = 1/(2*sigma^2);

m = size(Data,1);

index = StatePart.getElementPartition(currentState).index;
Input = Data(index,3);

for i =1:m
    %tempFuncOut = @(y) GaussianKernel(Data(i,1),y,gamma,[],'Normal');
    %Phi{i} = tempFuncOut;
    tempFuncOut = @(x,u) exp(1)*GaussianKernel(Data(i,2),x,gamma,[],'Normal')*GaussianKernel(Data(i,3),u,2*sigma^2*gamma,[],'Normal');
    Psi{i} = tempFuncOut;
end

%PhiValue = evalPhi(Phi,currentState*ones(m,1));

out = evalPsi(Psi,currentState*ones(m,1),Input*ones(m,1));


end


function PsiValue = evalPsi(Psi,x,u)

N = length(Psi);

if length(x) ~= N || length(u) ~= N
    error('Dimensions of Phi and x must agree');
end

PsiValue = zeros(N,1);

for i=1:N
    PsiValue(i) = Psi{i}(x(i),u(i));
end

end