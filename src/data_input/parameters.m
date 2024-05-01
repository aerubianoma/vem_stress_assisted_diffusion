function [lambda,mu,m0,m1,phi0,phi1,theta] = parameters()

    lambda = 1.e3;
    mu = 2*1.e2; % Multiply by 2 since the code is made for \mu = \mu/2.
    m0 = 1.e-1;
    m1 = 1.e-4;
    phi0 = 1;
    phi1 = 1;
    theta = 1.e-3;