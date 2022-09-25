function out = GaussianKernel(x,y,gamma,alpha,mode)


switch mode
    case 'Normal'
        if size(x) ~= size(y)
            error('The first two inputs must have the same dimension');
        end
        
        SumKernel = 0;
        
        if size(x,2) > 1
            L = size(x,2);
            for i =1:L
                for j =1:L
                    SumKernel = SumKernel + alpha(i)*alpha(j)*exp(-gamma*norm(x(:,i)-y(:,j))^2);
                end
            end
            
            out = sqrt(SumKernel);
        else
            out = exp(-gamma*norm(x-y)^2);
        end
    case 'MatrixFac'
        L = size(x,2);
        A = zeros(L,L);
        
        for i=1:L
            for j=1:L
                A(i,j) = exp(-gamma*norm(x(:,i)-y(:,j))^2);
            end
        end
        
        out = chol(A + 0.1*eye(L));
    otherwise
        error('This mode has not been implemented');
end



end

