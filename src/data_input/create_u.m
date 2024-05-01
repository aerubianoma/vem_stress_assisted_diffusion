function [u1,u2] = create_u()

    syms x y;
    
    %% For convergence analysis in L domain [-1,-1]^2\(0,1]^2
    
    %u1 = ((x^2+y^2)^(-1/6))*y;
    %u2 = ((x^2+y^2)^(-1/6))*(-x);

    %% For convergence analysis in square domain

    u1 = 0.2*cos(x)*sin(y)*x+0.2*x^2;
    u2 = 0.2*sin(x)*cos(y)*x+0.2*y^2;