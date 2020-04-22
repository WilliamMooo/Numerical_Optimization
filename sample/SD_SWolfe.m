%% ************************************************************
%  filename: SD_SWolfe
% 
%% *************************************************************
%% The steepest descent method + Strong Wolfe line search for solving
%
%  min f(x):= 0.5 ||Ax-b||^2
%

function [xsol,fsol,iter,ttime] = SD_SWolfe(x,OPTIONS,A,b)

if isfield(OPTIONS,'tol');       tol       = OPTIONS.tol;        end
if isfield(OPTIONS,'printyes');  printyes  = OPTIONS.printyes;   end
if isfield(OPTIONS,'maxiter');   maxiter   = OPTIONS.maxiter;    end

LSOPTIONS.printyes = 0;

LSOPTIONS.maxiter = 100;

lsopts.c1 = 1e-4;

lsopts.c2 = 0.9; 

lsopts.c2 = 0.9; 

alpha_max = 100;

normb = norm(b);

if  (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n \t   Steepest descent method for solving unconstrained minimization');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter    lstep      rnormg        obj     cratio    time');
end


tstart = clock;

xk = x;

Axkb = A*xk - b;

fk = 0.5*norm(Axkb)^2;

gk = A'*Axkb;

rnormg = norm(gk)/max(1,normb);

iter = 1; 

while (rnormg > tol) || (iter<=maxiter)
    
    %----------- calculate the direction -----------
    
    dk = -gk;
    
    cratio = gk'*dk;     % the curvature ratio
    
    alphak = lsearch_SWolfe(xk,fk,gk,dk,LSOPTIONS,lsopts,alpha_max,10,A,b);
    
    ttime = etime(clock,tstart);
     
    if (printyes)
        
        fprintf('\n %3.0d    %3.2e     %3.2e    %3.2e   %3.2e    %3.2f',iter,alphak,rnormg,fk,abs(cratio),ttime);
    end
    
    xnew = xk+alphak*dk;
    
    Axnewb = A*xnew - b;
    
    fnew = 0.5*norm(Axnewb)^2;
    
    gnew = A'*Axnewb;
          
    rnormg = norm(gnew)/max(1,normb);
       
    xk = xnew; fk = fnew;
    
    gk = gnew;
    
    iter = iter + 1;
end

xsol = xk;

fsol = fnew;

end