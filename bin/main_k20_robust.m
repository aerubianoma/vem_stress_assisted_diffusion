clc; 
close all;
clear variables;

format short g

% run time check
tic;

%% Mesh size

mesh = [8;12;16;20];

maxIt = size(mesh,1);

%% Neumann condition

bdNeumann = 'x==1 | y==1'; 

%% For error purposes

h = zeros(maxIt,1);
N_psp1 = zeros(maxIt,1);
N_psp2 = zeros(maxIt,1);

Err_u = zeros(maxIt,1);
Err_p = zeros(maxIt,1);
Err_zeta = zeros(maxIt,1);
Err_phi = zeros(maxIt,1);

Err_u_ratio = zeros(maxIt,1);
Err_p_ratio = zeros(maxIt,1);
Err_zeta_ratio = zeros(maxIt,1);
Err_phi_ratio = zeros(maxIt,1);

%% Physical parameters

[lambda,mu,m0,m1,phi0,phi1,theta] = parameters();

%% For $\lambda$

variation = [5.e1;1.e2;5.e2;1.e3;5.e3];

%% For $\mu$

%variation = 2*[5.e1;1.e2;5.e2;1.e3;5.e3];

%% For $\theta$

%variation = [5.e-4;1.e-4;5.e-3;1.e-3;5.e-2];

%% Save errors and $M$ constant

errors_V1_Q1 = zeros(maxIt,size(variation,2));
errors_V2_Q2 = zeros(maxIt,size(variation,2));
M_values = zeros(size(variation));

for c = 1:size(variation,1)

    %% Variable to check robustness change name corresponding to the parameter
    
    lambda = variation(c);
    variable_string = '$\lambda$';
    
    fprintf('constant value: %d\n', variation(c));
    
    %% Data for vem
    
    data_stress = mixed_concentration_active_stress_data_MBC(mu,lambda,phi0,phi1);
    data_diffusion = mixed_stress_assisted_diffusion_data_MBC(mu,lambda,theta,phi0,phi1,m0,m1);
    
    %% Obtain M parameter using optimization routine
    
    xM = fmincon(@(p) -norm(data_diffusion.M_1(p),'fro'),[0,0],[],[],[],[],[0,0],[1,1]);
    M = data_diffusion.M_1(xM);
    M = abs(M(1,1));
    
    %% Save M value
    
    M_values(c) = M;
    
    %% Stop tolerance for Piccard iteration
    
    tol = 1.e-8;
    
    %% Convergence rate experiment
    
    for k = 1:maxIt
        
        %% unit-square with different discretisations: 
    
        % nonconvex, polygonal, square, distortionpolygonal, kangaroo
        
        %% Load mesh and domain
    
        mesh_type = 'nonconvex';
        mesh_file = strcat(mesh_type, num2str(mesh(k)), '.mat');
        load(mesh_file);
        
        %% print actual mesh and # elements
    
        fprintf(strcat('mesh: ',mesh_file, '\n'));
        fprintf('Number of elements: %d\n', size(elem,1))
    
        %% Set boundary information
    
        bdStruct = setboundary(node,elem,bdNeumann);
    
        %% initial conditions:
        %% Coefficients of a polynom in \mathcal{M}_{1}
    
        info_diffusion = 1;
    
        %% Error for fixed point iteration
    
        Err_u = [];
        Err_p = [];
        Err_zeta = [];
        Err_phi = [];
        Err = [];
    
        it = [];
        it_ct = 0;
    
        %% Fixed point iteration
    
        while true
    
            %% Solve using vem
        
            %% k_1 = 2 for elasticity mixed formulation (mixed boundary conditions)
            
            [uh,ph,info_stress_next] = vem_mixed_concentration_active_stress_k0_MBC(node,elem,data_stress,bdStruct,phi0,phi1,info_diffusion);
            
            %% k_2 = 0 for transport mixed formulation (mixed boundary conditions) 
            
            [zetah,phih,info_diffusion_next] = vem_mixed_stress_assisted_diffusion_k0_MBC(node,elem,data_diffusion,bdStruct,m0,m1,info_stress_next);
        
            %% Record errors, all calculated in the scalled norms \mathcal{V}_1,Q_1,\mathcal{V}_2,Q_2
        
            % \mathcal{V}_1 error and norm
    
            [Err_u_it,u_it] = error_V1(node,elem,uh,info_stress_next,data_stress);
    
            % Q_1 error and norm
    
            [Err_p_it,p_it] = error_Q1(node,elem,ph,data_stress);
            
            % \mathcal{V}_2 and Q_2 error and norm
            
            [Err_zeta_it,Err_phi_it,zeta_it,phi_it] = error_V2_Q2_k0(node,elem,zetah,phih,info_diffusion_next,data_diffusion,M);
    
            %% Store errors for tables and graphs
    
            Err_u = [Err_u, Err_u_it];
            Err_p = [Err_p, Err_p_it];
            Err_zeta = [Err_zeta, Err_zeta_it];
            Err_phi = [Err_phi, Err_phi_it];
    
            %% Track iterations
    
            it_ct = it_ct+1;
            it = [it, it_ct];
    
            %% Stop criteria
    
            if size(it,2) > 1
                
                %% \norm{sol_{l} - sol_{l-1}} < tol for stop criteria
    
                [Err_u_stop,u_stop] = error_V1(node,elem,uh-uh_prev,info_stress_next,data_stress);
                [Err_L2_stop,p_stop] = error_Q1(node,elem,ph-ph_prev,data_stress);     
                [Err_zeta_stop,Err_phi_stop,zeta_stop,phi_stop] = error_V2_Q2_k0(node,elem,zetah-zetah_prev,phih-phih_prev,info_diffusion_next,data_diffusion,M);
    
                if (sqrt(u_stop.^2+p_stop.^2)+sqrt(zeta_stop.^2+phi_stop.^2)) < tol
    
                    fprintf('Total iterations: %d\n', it(end))
                    break
                
                end
            
            end
    
            %% Update solutions
    
            info_stress = info_stress_next;
            info_diffusion = info_diffusion_next;
    
            uh_prev = uh;
            ph_prev = ph;
            zetah_prev = zetah;
            phih_prev = phih;        
            
        end
    
        % For error plots and tables
    
        N_psp1(k) = length(uh);
        N_psp2(k) = length(zetah);  
        h(k) = 1/sqrt(size(elem,1));
    
        Err_zeta_ratio(k) = Err_zeta(end);
        Err_phi_ratio(k) = Err_phi(end);
        Err_u_ratio(k) = Err_u(end);
        Err_p_ratio(k) = Err_p(end);
    
    end

errors_V1_Q1(:,c) = sqrt(Err_u_ratio.^2+Err_p_ratio.^2);
errors_V2_Q2(:,c) = sqrt(Err_zeta_ratio.^2+Err_phi_ratio.^2);
    
end

total_Err = sqrt(errors_V1_Q1.^2+errors_V2_Q2.^2);

figure(1);

format short g

showrate_robust(h,total_Err,variation,M_values,variable_string,'$M$',1)

filename = strcat('../outputs/robust/',mesh_type,'/k20_',variable_string);
saveas(gcf, fullfile(strcat(filename, '.fig')));

toc