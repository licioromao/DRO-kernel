function out = GaussianKernel(x,y,gamma,alpha)

if size(x) ~= size(y)
    error('The first two inputs must have the same dimension');
end

SumKernel = 0;

if size(x,2) > 1
    L = size(x,2);
    for i =1:L
        for j =1:L
            if i ~= L && j ~= L
                SumKernel = SumKernel + alpha(i)*alpha(j)*exp(-gamma*norm(x(:,i)-y(:,j))^2);
            else
                if i == L && j ~=L
                    SumKernel = SumKernel + alpha(i)*alpha(j)*exp(-gamma*norm([210;210;1300]-y(:,j))^2);  
                end
                
                if i == L && j == L
                    SumKernel = SumKernel + alpha(i)*alpha(j);
                end
                
                if i ~= L && j ==L
                    SumKernel = SumKernel + alpha(i)*alpha(j)*exp(-gamma*norm(x(:,i) - [210;210;1300])^2);
                end
            end
        end
    end
    
    out = sqrt(SumKernel);
else
    out = exp(-gamma*norm(x-y)^2);
end


end

