%% *************************************************************
%  filename: test_Conjugate
%
%%  ****************** generate the problem *******************

clear;

addpath(genpath(pwd));

n =3000;

%% **************** to fix the random seed **********************
randstate = 10
randn('state',double(randstate));
rand('state',double(randstate)); 

%% *************************************************************
         
B = randn(n,n);

%% ************ one way to generate A **************************

A = B*B'+ 0.01*eye(n);

b = randn(n,1);

c = randn(1,1);

%% ************ Parameters for steepest descent method ***********

OPTIONS.tol = 1.0e-5;

OPTIONS.maxiter = 30000;

OPTIONS.printyes = 0;
 
x0 = zeros(n,1);       %% the starting point must be chosen as this !!!!

 tic
 res=b;
 [xsol,iter] = CG_scale (A,res,b,1e-5,3000);
 toc
