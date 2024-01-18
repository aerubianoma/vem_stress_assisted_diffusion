% Convergence test in 2D 
% Steady coupled stress-assisted diffusion in double perturbed saddle-point formulation 
% 
% Nonlinearities treated via (manual) fixed-point iterations
% 
% 
% Strong primal form: 
%                          - div(sigma) = f
% alpha*phi - div(M(sigma) * grad(phi)) = g
% 
% with 
% 
% sigma = 2*mu*eps(u) + lambda*div(u)*I - ell(phi)*I
% M     = m0*(1+exp(-m1*tr(sigma)))*I
% 
% CHECK BOUNDS underline(M) and overline(M)!!!
% 
% 
% and non-homogeneous mixed boundary conditions 
% 
% u = uD and phi = phiD on GammaD
% sigma*n = tN and M(sigma)*grad(phi).n = sN
% 
% The domain is the unit square, divided into GammaD = {x=0 U y=0} and GammaN = partialOmega\GammaD 
% 
% Both sub-problems are written in mixed form, using the new variables 
% 
% tp   = - lambda*div(u) + ell(phi)
% zeta = M(sigma)*grad(phi) 
% 
% Strong mixed form: 
%       - div(2*mu*eps(u)-tp*I) = f
%        -div(u) - 1/lambda*tp  = -ell(phi)
% M^{-1}(u,tp)*zeta - grad(phi) = 0
%         alpha*phi - div(zeta) = g                                           
% 
% 
% u  tp zeta phi are in 
% H1 L2 Hdiv L2
% 
% Unit square, manufactured solutions

function M1_l = M1_local(m0,m1,mu,c_w,c_r,hK,xK,yK)
    
    [M1_11,M1_12,M1_22] = compute_M1(m0,m1,mu,c_w,c_r,hK,xK,yK);

    function val = M1(p)
        x = p(:,1); 
        y = p(:,2);
        row1 = [M1_11(x,y) + 0*x, M1_12(x,y) + 0*x];
        row2 = [M1_12(x,y) + 0*x, M1_22(x,y) + 0*x];
        val = [row1; row2];
    end


    M1_l = struct('M1', @M1);    
end

    function [M1_11,M1_12,M1_22] = compute_M1(m0,m1,mu,c_w,c_r,hK,xK,yK)

    syms x y;

    %% Basis for $P_2(E)$ %%
    
    m1P2 = 1+0*x;                  
    m2P2 = (x-xK)./hK+0*x;             
    m3P2 = (y-yK)./hK+0*x;             
    m4P2 = (x-xK).^2/hK^2+0*x;         
    m5P2 = (x-xK).*(y-yK)./hK^2+0*x;   
    m6P2 = (y-yK).^2./hK^2+0*x;  
    
    %% Basis for $P_1(E)$

    m1P1 = 1+0*x;
    m2P1 = (x-xK)./hK; 
    m3P1 = (y-yK)./hK ; 

    % previous iteration solution

    w1 = c_w(1)*[m1P2,0*x]+c_w(2)*[m2P2,0*x]+c_w(3)*[m3P2,0*x]+c_w(4)*[m4P2,0*x]+c_w(5)*[m5P2,0*x]+c_w(6)*[m6P2,0*x];
    w2 = c_w(7)*[0*x,m1P2]+c_w(8)*[0*x,m2P2]+c_w(9)*[0*x,m3P2]+c_w(10)*[0*x,m4P2]+c_w(11)*[0*x,m5P2]+c_w(12)*[0*x,m6P2];

    w1 = w1(1,1);
    w2 = w2(1,2);
    
    r = c_r(1)*m1P1+c_r(2)*m2P1+c_r(3)*m3P1;

    % derivative

    w1x = diff(w1,x);
    w2y = diff(w2,y);

    M_11 = m0*exp(-m1*(mu*(w1x+w2y)-2*r)) + 0*x;
    M_12 = 0 + 0*x;
    M_22 = m0*exp(-m1*(mu*(w1x+w2y)-2*r)) + 0*x;
    detM = M_11*M_22-M_12.^2;

    M1_11 = M_22./detM;
    M1_12 = -M_12./detM;
    M1_22 = M_11./detM;

    M1_11 = matlabFunction(M1_11,'Vars',{x,y});  
    M1_12 = matlabFunction(M1_12,'Vars',{x,y}); 
    M1_22 = matlabFunction(M1_22,'Vars',{x,y});

end