function l_l = l_local_k_1_battery(phi0,phi1,c_phi,hK,xK,yK)
    
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

    %% Hill's function approach

    %l = phi0 + phi^2/(phi1+phi^2);

    %% Taralovas approach

    l = phi0*phi;

    l = matlabFunction(l,'Vars',{x,y});  

end