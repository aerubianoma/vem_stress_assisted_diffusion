function pde = mixed_concentration_active_stress_data_MBC_battery(mu,lambda)
% PDE: u = [u1,u2]
% coupled term \varphi
%     -mu*\vdid(\beps(u) - pI) = f,  \in \Omega
%     p = -\lambda*\vdiv u + \ell(\varphi),  \in \Omega
%     u = g_D, in \Gamma_D
%     (\beps(u) - pI) \cdot n = g_N in \Gamma_N 

    % load source
    function val = f(p)
        x = p(:,1); y = p(:,2);
        val = [0*x,0*x]; % Taralova's thesis % % 75.*[x+0*x,y+0*x]; % Ricardo's paper %
    end

    % clamped
    function val = u_D(p)
        x = p(:,1); y = p(:,2);
        val = [0+0*x,0+0*x];
    end

    function val = sigma_N(p)
        x = p(:,1); y = p(:,2);
        val = 0 + 0*x; %0 + 0*x;
    end

pde = struct('f',@f,'u_D',@u_D, 'sigma_N',@sigma_N, 'mu',mu, 'lambda',lambda);

end