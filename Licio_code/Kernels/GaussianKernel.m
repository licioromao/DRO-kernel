function out = GaussianKernel(x,y,gamma,alpha,mode,input_vector)



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
        
        out = chol(A + 0.01*eye(L));
    case 'KME'
        number_of_points = size(x,1);
        
        if number_of_points ~= size(y,1)
            error('The size of the first two inputs must coincide');
        end

        A = zeros(number_of_points,number_of_points);
        
        for i=1:number_of_points
            state_1 = x(i,:);
            input_1 = input_vector(i,:);

            for j=1:number_of_points
                state_2 = x(j,:);
                input_2 = input_vector(j,:);

                A(i,j) = exp(-gamma*norm(state_1 - state_2)^2) + ...
                            exp(-gamma*norm(input_1 - input_2)^2);
            end
        end
        
        out = A;

    case 'KME_2'
        state_1 = x;
        state_2 = y;
        input_1 = input_vector(1,1);
        input_2 = input_vector(1,2);

        out =  exp(1)*exp(-gamma*norm(state_1 - state_2)^2) + ...
            0.001*linear_splines(input_1,input_2);
    otherwise
        error('This mode has not been implemented');
end

    
end

function out = linear_splines(input_1,input_2)

    out = 1 + input_1*input_2 + input_1*input_2*min(input_1,input_2)...
            - ((input_1+input_2)/2)*min(input_1,input_2)^2 + (1/3)*min(input_1,input_2)^3;

end