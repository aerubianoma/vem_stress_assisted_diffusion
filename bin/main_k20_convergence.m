clc; 
close all;
clear variables;

%% Mesh size

% For kangaroo

mesh = [2;3;4;5];

% Other cases

%mesh = [8;12;16;20];

maxIt = size(mesh,1);

%% Neumann condition

bdNeumann = 'x==1 | y==1'; 

%% For error purposes

h = zeros(maxIt,1);
N_psp1 = zeros(maxIt,1);
N_psp2 = zeros(maxIt,1);
total_it = zeros(maxIt,1);
time_it = zeros(maxIt,1);

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

%% Data for vem

data_stress = mixed_concentration_active_stress_data_MBC(mu,lambda,phi0,phi1);
data_diffusion = mixed_stress_assisted_diffusion_data_MBC(mu,lambda,theta,phi0,phi1,m0,m1);

%% Obtain M parameter using optimization routine

xM = fmincon(@(p) -norm(data_diffusion.M_1(p),'fro'),[0,0],[],[],[],[],[0,0],[1,1]);
M = data_diffusion.M_1(xM);
M = abs(M(1,1));

fprintf('Value of M: %d\n', M)

%% Stop tolerance for Piccard iteration

tol = 1.e-8;

%% Convergence rate experiment

for actual_mesh = 1:maxIt
    
    tic;
    
    %% unit-square with different discretisations: 

    % nonconvex, polygonal, square, distortionpolygonal, crossed, kangaroo
    
    % CHANGE THIS FOR EXAMPLE 1 TABLE/GRAPH GENERATION! (line 71)

    %% Load mesh and domain

    mesh_type = 'kangaroo';
    mesh_file = strcat(mesh_type, num2str(mesh(actual_mesh)), '.mat');
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

    N_psp1(actual_mesh) = length(uh)+length(ph);
    N_psp2(actual_mesh) = length(zetah)+length(phih);  
    h(actual_mesh) = 1/sqrt(size(elem,1));
    total_it(actual_mesh) = it(end);
    time_it(actual_mesh) = toc;

    Err_zeta_ratio(actual_mesh) = Err_zeta(end);
    Err_phi_ratio(actual_mesh) = Err_phi(end);
    Err_u_ratio(actual_mesh) = Err_u(end);
    Err_p_ratio(actual_mesh) = Err_p(end);

    %% Projection for plotting
    
    [zetahI,phihI,nodeI_zetaphi,elemI_zetaphi] = projection_V2_k0(node,elem,zetah,phih,info_diffusion_next,uh*0);
    [phI,nodeI_p,elemI_p] = projection_Q1(node,elem,ph,1,uh*0);

    % Plots
    
    figure(1)
    showresult_u_p(node,elem,nodeI_p,elemI_p,data_stress.uexact,uh,data_stress.pexact,phI);
    stress_file = strcat('stress_',num2str(mesh(actual_mesh)),'.fig');
    saveas(gcf, fullfile('../outputs/graphs/',mesh_type,'/k20', stress_file));

    figure(2); 
    showresult_zeta_phi(nodeI_zetaphi,elemI_zetaphi,data_diffusion.zeta_exact,zetahI,data_diffusion.phi_exact,phihI);  
    diffusion_file = strcat('diffusion_',num2str(mesh(actual_mesh)),'.fig');
    saveas(gcf, fullfile('../outputs/graphs/',mesh_type,'/k20', diffusion_file));
    
    drawnow;   
    
end

Err_V1_Q1 = sqrt(Err_u_ratio.^2+Err_p_ratio.^2);
Err_V2_Q2 = sqrt(Err_zeta_ratio.^2+Err_phi_ratio.^2);

total_Err = sqrt(Err_V1_Q1.^2+Err_V2_Q2.^2);

filename = strcat('../outputs/tables/',mesh_type,'/k20_');

figure(5);
showrateh(h,total_Err,[0 0.4470 0.7410],'$||(\mathbf{u}-\mathbf{u}_h,\tilde{p}-\tilde{p}_h,\mathbf{\zeta}-\mathbf{\zeta}_h,\varphi-\varphi_h)||_{\mathbf{V}_1\times Q_1\times \mathbf{V}_2\times Q_2}$')
saveas(gcf, fullfile(strcat(filename, 'rate_graph.fig')));

total_rates = zeros(size(total_Err));

for i=1:size(total_rates,1)-1

    total_rates(i+1) = log(total_Err(i+1)/total_Err(i))/log(h(i+1)/h(i));

end

colname = {'#Dof', 'h', 'e', 'r', 'It', 't'};
disptable(strcat(filename,'table.txt'),colname,N_psp1+N_psp2,[],h,'%0.2e',total_Err,'%0.2e',total_rates,'%1.2e', total_it, [], time_it, '%1.2e');

disp(strcat('Total run time: ',num2str(sum(time_it))))