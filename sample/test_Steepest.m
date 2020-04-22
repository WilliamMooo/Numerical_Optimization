%% *************************************************************
%  filename: test_Steepest
%
%%  ****************** generate the problem *******************
addpath(genpath('solvers'));

n = 100;

B = randn(n,n);

A = (B+B')/2;

root = eig(A);

lamda = root(1);

if(lamda<0)
    A = A - 2*lamda*eye(n);
end

b = randn(n,1);

%% ************ Parameters for steepest descent method *******

OPTIONS.tol = 1.0e-6;

OPTIONS.maxiter = 5;

OPTIONS.printyes = 1;
 
x0 = ones(n,1);   % starting point

[xsol,fsol,iter,time] = SD_SWolfe(x0,OPTIONS,A,b);

% [xsol,fsol,iter,time] = GS_method(x0,OPTIONS,A,b,c);
% 
% [xsol,fsol,iter,time] = sGS_method(x0,OPTIONS,A,b,c);

%[xsol,fsol,iter,time] = Steepest_exact(x0,OPTIONS,A,b,c);
