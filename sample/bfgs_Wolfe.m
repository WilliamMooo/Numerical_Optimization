%% *************************************************************
%  filename: bfgs_Wolfe
%% *************************************************************
%% The bfgs+Wolfe code for solving the following unconstrained minimization
%  min f(x)
%  where f is continuously differentiable
%%

function [xsol,fsol,iter,ttime] = bfgs_Wolfe(x,OPTIONS)

if isfield(OPTIONS,'tol');       tol       = OPTIONS.tol;        end
if isfield(OPTIONS,'printyes');  printyes  = OPTIONS.printyes;   end
if isfield(OPTIONS,'maxiter');   maxiter   = OPTIONS.maxiter;    end

LSOPTIONS.printyes = 0;

LSOPTIONS.maxiter = 100;

lsopts.c1 = 1e-4;

lsopts.c2 = 0.9; 

alpha_max =100;

eta = 1.0e-8;

n = size(x,1);

if  (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n \t   the bfgs method for solving unconstrained minimization');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter    fobj         normg        lstep          sy       time ');
end

%% ********** to yield an initial inverse Hessian ****************

obj_list = zeros(maxiter,1);

iter = 1;

tstart = clock;

[fobj,g] = objfun(x);  

dir = -g;

alpha = lsearch_SWolfe(x,fobj,g,dir,LSOPTIONS,lsopts,alpha_max,50);

xnew = x + alpha*dir;

[fnew,gnew] = objfun(xnew);  

dx = xnew - x;  dg = gnew - g;  

g = gnew;   x = xnew;  fobj = fnew;

dxg = dx'*dg;

scale = dxg/norm(dg)^2;    % H0 = scale*eye(n)

Hessian = scale*eye(n);

dir = -Hessian*g;

normg = norm(g);

obj_list(iter)=fobj;

%% ******************* Main Loop ********************************

while (normg>tol && iter <= maxiter)

    alpha = lsearch_SWolfe(x,fobj,g,dir,LSOPTIONS,lsopts,alpha_max,100);

    xnew = x + alpha*dir;
    
    ttime = etime(clock,tstart);
    
    [fnew,gnew] = objfun(xnew);
 
    dx = xnew - x;  dg = gnew - g;
    
    g = gnew;   x = xnew;  fobj = fnew;
    
    dxg = dx'*dg;  normg = norm(g);
    
    if (printyes)
        
        fprintf('\n %3.0d    %3.2e      %3.2e      %3.2e     %3.2e    %3.2f',iter,fobj,normg,alpha,dxg,ttime);
        
    end
    
    %% ***************** to yield a new direction *********************
    
    dxnorm = norm(dx);   dgnorm = norm(dg);
     
    if (dxg <= eta*dxnorm*dgnorm)
  
         dxg = dx'*dg;
        
        scale = dxg/norm(dg)^2;    % H0 = scale*eye(n)
        
        dir = -scale*g;
        
    else     
        
        w = Hessian*dg;
        
        zhz = dg'*w;
        
        c2 = 1.0/dxg;
    
        c1 = (1.0+c2*zhz)*c2;
        
        Hessian = Hessian -c2*(dx*w'+w*dx')+c1*(dx*dx');
        
        dir = -Hessian*g;
   
    end
    
    iter = iter + 1;
    
    obj_list(iter) = fobj;
    
    if (iter>=500)&& abs(obj_list(iter)-obj_list(iter-10))<=1e-2*tol
        
        xsol = x;  fsol=fobj;
        
        return;
    end

end

xsol = x;

fsol = fobj;

end

