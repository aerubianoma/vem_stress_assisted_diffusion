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
% sigma = _phi()2*mu*eps(u) + lambda*div(u)*I - ell(phi)*I
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

function l_l = l_local_k_1(phi0,phi1,c_phi,hK,xK,yK)
    
    l = compute_l(phi0,phi1,c_phi,hK,xK,yK);

    function val = g(p)
        x = p(:,1); y = p(:,2);
        val = l(x,y)+0*x;
    end

    l_l = struct('g',@g);
 
end

    function l = compute_l(phi0,phi1,c_phi,hK,xK,yK)

    syms x y;

    m1P1 = 1+0*x;
    m2P1 = (x-xK)./hK; 
    m3P1 = (y-yK)./hK ; 

    % previous iteration solution
    
    phi = c_phi(1)*m1P1+c_phi(2)*m2P1+c_phi(3)*m3P1;

    l = phi0 + phi^2/(phi1+phi^2);

    l = matlabFunction(l,'Vars',{x,y});  

end