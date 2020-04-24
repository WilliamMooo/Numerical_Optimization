%% ***************************************************************
%  filename: lsearch_SWolfe
% 
%% ***************************************************************
%% To search alpha* in (0,alpha_max) to satisfy the strong Wolfe 

function [alpha_star,iter] = lsearch_SWolfe(xk,fk,gk,dk,OPTIONS,lsopts,alpha_max,lsmax)

if isfield(OPTIONS,'printyes');  printyes  = OPTIONS.printyes;   end
if isfield(OPTIONS,'maxiter');   lsmax     = OPTIONS.maxiter;    end

c1 = lsopts.c1;

c2 = lsopts.c2;

cratio = gk'*dk;

if (printyes)
    fprintf('\n *****************************************************');
    fprintf('******************************************');
    fprintf('\n \t  line search with strong Wolfe-rule ');
    fprintf('\n ****************************************************');
    fprintf('*******************************************');
    fprintf('\n  iter    alpha_old       alpha      alpha_max');
end

alpha_old = 0;

fold = fk;

alpha = alpha_max/2;

%% ********************* Main Loop ********************************

for i=1:lsmax
    
    if (printyes)
        
        fprintf('\n %3.0d    %3.2e     %3.2e     %3.2e',i,alpha_old,alpha,alpha_max);
        
    end
    
    xnew = xk+alpha*dk;
    
    [fnew,gnew] = objfun(xnew);
    
    if (fnew> fk + c1*alpha*cratio) || (i>1 && fnew>=fold)
        
        alpha_star = zoom(alpha_old,alpha,fold,xk,fk,dk,cratio,lsopts);
        
        iter = i;
        
        return;
    end
    
    cratio_new = gnew'*dk;
    
    if abs(cratio_new)<=-c2*cratio  %% successful!!!
        
        alpha_star = alpha;
        
        iter = i;
        
        return;
        
    elseif cratio_new>=0
        
        alpha_star = zoom(alpha,alpha_old,fnew,xk,fk,dk,cratio,lsopts);
        
        iter = i;
        
        return;
    end
    
    alpha_old = alpha;  fold = fnew;
            
    alpha = alpha_old+(alpha_max-alpha_old)/2;
           
end

alpha_star = alpha;

iter = i;

%% ********************* zoom function *****************************
%%
function alpha_star = zoom(lalpha,halpha,flow,xk,fk,dk,cratio,lsopts)

c1 = lsopts.c1;

c2 = lsopts.c2;

flag = 1;

while (flag)
    
    talpha = lalpha+(halpha-lalpha)/2;
    
    xt = xk + talpha*dk;
    
    [ft,gt] = objfun(xt);
    
    if (ft> fk + c1*talpha*cratio) || (ft>=flow)
        
        halpha = talpha;
        
    else
        tcratio = gt'*dk;
        
        if abs(tcratio)<=-c2*cratio
            
            alpha_star = talpha;
            
            return;
        end
        
        if tcratio*(halpha-lalpha)>=0
            
            halpha = lalpha;
            
            lalpha = talpha;  
            
            flow = ft;
        end
    end
end


%% ******************* Objection function *****************************
%%

% function [fx,gradf] = objfun(x,A,b)
% 
% Axb = A*x-b;
% 
% fx = 0.5*norm(Axb)^2;
% 
% if nargout>1
%    
%   gradf = A'*Axb;
%   
% end