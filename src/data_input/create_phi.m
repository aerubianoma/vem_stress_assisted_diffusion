function phi = create_phi()

    syms x y;
    
    %% For convergence analysis in L domain [-1,1]^2\(0,1]^2

    %phi = sin(pi*x)+sin(pi*y)*exp(x+y);

    %% For convergence analyisis in square domain

    phi = x^2+y^2+sin(pi*x)+cos(pi*y);