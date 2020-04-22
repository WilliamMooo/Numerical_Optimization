%% **********************************************************************
%  filename: CG_scale
%
%% **********************************************************************
%% conjugate gradient method for solving the system of linear equations
% 
%   Hx = b  with  H = U + AVA'
%
%  x0£ºthe starting point
%
%  res£ºthe residual of Hx = b at x0, i.e., res = b-H(x0) 
%%
%% **********************************************************************************  

function [x,iter,solve_ok] = CG_scale(A,res,b,tol,maxit)
 
m = size(b,1); resnrm =[];

solve_ok = 1;
 
tiny = 1.0e-16;
 
stagnate_check = 20;

d = zeros(m,1);   Ad = zeros(m,1);

%% *************** Initialization part **************************

x =  zeros(m,1);       % Such a starting point is crucial !!!

r = res;
 
err = norm(r);
 
resnrm(1) = err;  minres = err;
 
z = r; 
 
rho_old = r'*z;
 
tau_old = norm(z);
 
theta_old = 0;

%% ******************* Main Loop ********************************

for iter = 1:maxit

    Az = A*z; 
    
    sigma = z'*Az;
        
    if (abs(sigma)<tiny)   %% in this case z=0 since A is positive definite
        
        solve_ok = -1;
        
        break;
        
    else
        
        alfa = rho_old/sigma;
        
        r = r - alfa * Az; 
        
    end
    
    u = r;
    
    theta = norm(u)/tau_old;
    
    c=1/sqrt(1+theta^2);
    
    tau = tau_old*theta*c;
    
    gam = (c*theta_old)^2;
    
    eta = alfa/(1+theta^2);
    
    d =  gam*d + eta*z;
    
    x = x + d;
    
    %%-------------------- stopping conditions  ---------------------------
    
    Ad = gam*Ad + eta*Az;
    
    res = res - Ad;
    
    err = norm(res);
    
    resnrm(iter+1) = err;
    
    [~,indictor] = find(err < minres);
    
    minres(indictor) = err(indictor);
    
    if (err< tol); 
       
        break; 
    end
    
    if (iter > stagnate_check) && (iter > 10)
        
        ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter);
        
        if (min(ratio) > 0.997) && (max(ratio) < 1.003)
            
            solve_ok = -2;
            
            break;
        end
    end
    
    %%------------------------------------------------------------------
    
    rho = r'*u;
    
    beta = rho/rho_old;
    
    z = u + beta*z;
    
    rho_old = rho;
    
    tau_old = tau;
    
    theta_old = theta;
end
 
if (iter == maxit); solve_ok = -3; end
 
end    

%%%%%%%%%%%%%%%%%%%%%% End of conjugate_gradient.m  %%%%%%%%%%%%%%%%%%%%%%%   