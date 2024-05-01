function pde = mixed_concentration_active_stress_data_MBC(mu,lambda,phi0,phi1)
% PDE: u = [u1,u2]
% coupled term \varphi
%     -mu*\vdid(\beps(u) - pI) = f,  \in \Omega
%     p = -\lambda*\vdiv u + \ell(\varphi),  \in \Omega
%     u = g_D, in \Gamma_D
%     (\beps(u) - pI) \cdot n = g_N in \Gamma_N 

% --------- given by the symbolic computation ------
[u1,u2,pe,f1,f2,l,u1x,u1y,u2x,u2y] = compute_rhs(mu,lambda,phi0,phi1);

mu = mu;
lambda = lambda;

% exact solution
    function val = uexact(p)
        x = p(:,1); y = p(:,2);
        val = [u1(x,y)+0*x,u2(x,y)+0*x];
    end
    function val = pexact(p)
        x = p(:,1); y = p(:,2);
        val = pe(x,y)+0*x;
    end
% right side hand function
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [f1(x,y)+0*x,f2(x,y)+0*x];
    end
    function val = g(p)
        x = p(:,1); y = p(:,2);
        val = l(x,y)+0*x;
    end

% Dirichlet
    function val = g_D(p)
        val = uexact(p);
    end

% Neumann boundary conditions ( right side hand )
    function val = g_N(p)
        val = mu*Eu(p)-pexact(p)*eye(2);
    end

% Du
    function val = Du(p)
        x = p(:,1); y = p(:,2);
        Du1 = [u1x(x,y)+0*x, u1y(x,y)+0*x];
        Du2 = [u2x(x,y)+0*x, u2y(x,y)+0*x];
        val = [Du1, Du2];
    end
%Eu
function val = Eu(p)
        x = p(:,1); y = p(:,2);
        Eu1 = [u1x(x,y)+0*x, 0.5*(u1y(x,y)+u2x(x,y))+0*x];
        Eu2 = [0.5*(u1y(x,y)+u2x(x,y))+0*x, u2y(x,y)+0*x];
        val = [Eu1, Eu2];
end

pde = struct('uexact',@uexact, 'pexact',@pexact, 'f',@f, 'g',@g, 'g_D',@g_D, 'Du',@Du, 'Eu',@Eu, 'g_N',@g_N, 'mu',mu, 'lambda',lambda);
end

function [u1,u2,pe,f1,f2,l,u1x,u1y,u2x,u2y] = compute_rhs(mu,lambda,phi0,phi1)
    
    syms x y;

    %% manufacture solution
    
    phi = create_phi;

    [u1,u2] = create_u;    

    l = phi0 + phi.^2/(phi1+phi.^2);

    % derivative
    u1x = diff(u1,x);   u1y = diff(u1,y);
    u2x = diff(u2,x);   u2y = diff(u2,y);
    pe = -lambda*(u1x+u2y)+l;
    u1xx = diff(u1x,x); u1yy = diff(u1y,y); u1yx = diff(u1y,x);
    u2xx = diff(u2x,x); u2yy = diff(u2y,y); u2xy = diff(u2x,y);
    px = diff(pe,x);    py = diff(pe,y);
    
    % eps(u)
    epsu1 = u1xx + 0.5*(u1yy+u2xy);
    epsu2 = u2yy + 0.5*(u2xx+u1yx);

    f1 = -mu*epsu1 + px;
    f2 = -mu*epsu2 + py;

    % convert to anonymous functions
    
    u1 = matlabFunction(u1,'Vars',{x,y});
    u2 = matlabFunction(u2,'Vars',{x,y});
    pe = matlabFunction(pe,'Vars',{x,y});
    f1 = matlabFunction(f1,'Vars',{x,y});
    f2 = matlabFunction(f2,'Vars',{x,y});
    l = matlabFunction(l,'Vars',{x,y});
    u1x = matlabFunction(u1x,'Vars',{x,y});
    u2x = matlabFunction(u2x,'Vars',{x,y});
    u1y = matlabFunction(u1y,'Vars',{x,y});
    u2y = matlabFunction(u2y,'Vars',{x,y});
end