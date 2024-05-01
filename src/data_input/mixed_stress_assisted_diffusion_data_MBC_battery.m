function pde = mixed_stress_assisted_diffusion_data_MBC_battery(theta,mu)

    %% Strong form of the equation %%
    
    % \zeta = \mathcal{M}(w,r)*\nabla(\varphi)
    % \theta*\varphi - \vdiv \zeta = g
    % \varphi = \varphi_D \in \Gamma_D
    % \zeta \cdot n = \zeta_N \in \Gamma_N
    
    % net diffusive source
    
    function val = g(p)
        x = p(:,1); 
        y = p(:,2);
        val = 0 + 0*x;
    end

    % max concentration in \Gamma_D

    function val = phi_D(p)
        x = p(:,1); 
        y = p(:,2);
        val = 2.29e4 + 0*x; % 3.497e-6 + 0*x; % Ricardo's paper % % 3.497*rand*10^(-5) + 0*x; % random distribution of concentration %
    end

    function val = zeta_N(p)
        x = p(:,1); 
        y = p(:,2);
        val = -2 + 0*x; % 3.497e-6 + 0*x; % Ricardo's paper % % 3.497*rand*10^(-5) + 0*x; % random distribution of concentration %
    end


pde = struct( 'g', @g, 'phi_D', @phi_D, 'zeta_N', @zeta_N ,'theta', theta, 'mu', mu);

end
