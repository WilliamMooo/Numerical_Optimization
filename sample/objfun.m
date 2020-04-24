%% ***************************************************************
%  Filename: objfun
%
%% **************************************************************
function [fval,grad,H] = objfun(x)

n = size(x,1);

grad = zeros(n,1);

H = zeros(n,n);

fval = (x(2) - x(1)^2)^2 + (1 - x(1))^2;

for i = 2:floor(n/2)
    fval = fval + (x(2*i) - x(2*i-1)^2)^2 + (1 - x(2*i-1))^2;
end

if nargout==2
    
    fval = (x(2) - x(1)^2)^2 + (1 - x(1))^2;
    
    for i = 2:floor(n/2)
        
        fval = fval + (x(2*i) - x(2*i-1)^2)^2 + (1 - x(2*i-1))^2;
    end
    
    for i = 1:n
        if mod(i,2) ~= 0
            grad(i) = -2*(1 - x(i)) - 4*x(i)*(x(i+1) - (x(i))^2);
        else
            grad(i) = 2*(x(i) - (x(i-1))^2);
        end
    end
    
elseif nargout>2
    
    fval = (x(2) - x(1)^2)^2 + (1 - x(1))^2;
    
    for i = 2:floor(n/2)
        fval = fval + (x(2*i) - x(2*i-1)^2)^2 + (1 - x(2*i-1))^2;
    end
    for i = 1:n
        
        if mod(i,2) ~= 0
            grad(i) = -2*(1 - x(i)) - 4*x(i)*(x(i+1) - (x(i))^2);
            H(i,i) = 12*x(i)^2 - 4*x(i+1) + 2;
            H(i,i+1) = -4*x(i);
        else
            grad(i) = 2*(x(i) - (x(i-1))^2);
            H(i,i) = 2;
            H(i,i-1) = -4*x(i-1);
        end
    end
end
end