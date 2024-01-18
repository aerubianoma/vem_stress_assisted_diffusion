function pde = mixed_stress_assisted_diffusion_data_MBC(mu,lambda,theta,phi0,phi1,m0,m1)

    %% Strong form of the equation %%
    
    % \zeta = \mathcal{M}(w,r)*\nabla(\varphi)
    % \theta*\varphi - \vdiv \zeta = g
    % \varphi = \varphi_D \in \Gamma_D
    % \zeta \cdot n = \zeta_N \in \Gamma_N
    
    %M = [y^2+1, -x*y;
    %     -x*y , x^2+1]
    
    %% Given by the symbolic computation %%
    
    [phie,M_11,M_12,M_22,M1_11,M1_12,M1_22,Mgradphix,Mgradphiy,rhs,div_zeta] = compute_rhs(mu,lambda,theta,phi0,phi1,m0,m1);
    
    % coefficient \mathcal{M}(w,r)

    function val = M(p)
        x = p(:,1); 
        y = p(:,2);
        row1 = [M_11(x,y) + 0*x, M_12(x,y) + 0*x];
        row2 = [M_12(x,y) + 0*x, M_22(x,y) + 0*x];
        val = [row1; row2];
    end

   function val = M_1(p)
        x = p(:,1); 
        y = p(:,2);
        row1 = [M1_11(x,y) + 0*x, M1_12(x,y) + 0*x];
        row2 = [M1_12(x,y) + 0*x, M1_22(x,y) + 0*x];
        val = [row1; row2];
    end

    % exact solution
    
    function val = zeta_exact(p)
        x = p(:,1); 
        y = p(:,2);
        val = [Mgradphix(x,y) + 0*x,Mgradphiy(x,y) + 0*y];
    end
    
    function val = g(p)
        x = p(:,1); 
        y = p(:,2);
        val = rhs(x,y) + 0*x;
    end

    function val = phi_exact(p)
        x = p(:,1); 
        y = p(:,2);
        val = phie(x,y) + 0*x;
    end

    function val = div_z(p)
        x = p(:,1); 
        y = p(:,2);
        val = div_zeta(x,y) + 0*x;
    end

pde = struct('zeta_exact', @zeta_exact, 'g', @g, 'M', @M, 'M_1', @M_1, 'phi_exact', @phi_exact, 'theta', theta, 'mu', mu, 'div_z',@div_z);

end

function [phi,M_11,M_12,M_22,M1_11,M1_12,M1_22,Mgradphix,Mgradphiy,g,div_zeta] = compute_rhs(mu,lambda,theta,phi0,phi1,m0,m1)

    syms x y;

    % exact solution

    phi = create_phi;

    [u1,u2] = create_u;   

    l = phi0 + phi^2/(phi1+phi^2);
    
     % derivative

    u1x = diff(u1,x);
    u2y = diff(u2,y);

    phix = diff(phi,x);     
    phiy = diff(phi,y);

    p = -lambda*(u1x+u2y) + l;

    M_11 = m0*exp(-m1*(mu*(u1x+u2y)-2*p)) + 0*x;
    M_12 = 0 + 0*x;
    M_22 = m0*exp(-m1*(mu*(u1x+u2y)-2*p)) + 0*x;
    detM = M_11*M_22-M_12.^2;

    M1_11 = M_22./detM;
    M1_12 = -M_12./detM;
    M1_22 = M_11./detM;

    Mgradphix = M_11.*phix + M_12.*phiy;
    Mgradphiy = M_12.*phix + M_22.*phiy;

    % g = -\vdiv(M*\nabla(\varphi))+ \theta*\varphi

    g = -(diff(Mgradphix,x) + diff(Mgradphiy,y)) + theta*phi;

    div_zeta = diff(Mgradphix,x) + diff(Mgradphiy,y);

    % convert to anonymous functions
    
    phi = matlabFunction(phi,'Vars',{x,y});    

    g = matlabFunction(g,'Vars',{x,y});

    M_11 = matlabFunction(M_11,'Vars',{x,y});  
    M_12 = matlabFunction(M_12,'Vars',{x,y}); 
    M_22 = matlabFunction(M_22,'Vars',{x,y}); 

    M1_11 = matlabFunction(M1_11,'Vars',{x,y});  
    M1_12 = matlabFunction(M1_12,'Vars',{x,y}); 
    M1_22 = matlabFunction(M1_22,'Vars',{x,y});

    Mgradphix = matlabFunction(Mgradphix,'Vars',{x,y});  
    Mgradphiy = matlabFunction(Mgradphiy,'Vars',{x,y});
    div_zeta = matlabFunction(div_zeta,'Vars',{x,y});

end
