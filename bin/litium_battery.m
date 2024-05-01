clc; 
close all;
clear variables;

% run time check
tic;

%% Neumann condition

bdNeumann1 = '(x.^2+y.^2)==5.^2'; %'x==1 | y==1 | x==0 | y==0';%
bdNeumann2 = '(x.^2+y.^2)==1.^2'; % '(x.^2+y.^2)==1.^2'; %

%% Physical parameters

[lambda,mu,m0,m1,phi0,phi1,theta] = parameters_battery();

%% Constant stress-assisted coefficient

m0 = 1;
m1 = 0;

%% Data for vem

data_stress = mixed_concentration_active_stress_data_MBC_battery(mu,lambda);
data_diffusion = mixed_stress_assisted_diffusion_data_MBC_battery(theta,mu);

%% Obtain M parameter using optimization routine

M = 1;

fprintf('Value of M: %d\n', M)

%% Stop tolerance for Piccard iteration

tol = 1.e-8;
    
%% Load mesh and domain

mesh_file = 'meshes/circle_circle_10000.mat';
load(mesh_file);

%% print actual mesh and # elements

fprintf(strcat('mesh: ',mesh_file, '\n'));
h = size(elem,1);
fprintf('Number of elements: %d\n', h)

%% Set boundary information

bdStruct1 = setboundary_curve(1/sqrt(sqrt(h)),node,elem,bdNeumann1);
bdStruct2 = setboundary_curve(1/sqrt(sqrt(h)),node,elem,bdNeumann2);

showmesh(node,elem);

%% initial conditions:
%% Coefficients of a polynom in \mathcal{M}_{1}

info_diffusion = [1;0;0];

%% Tracking iterations

it = [];
it_ct = 0;

%% Fixed point iteration

while true

    %% Solve using vem

    %% k_1 = 2 for elasticity mixed formulation (mixed boundary conditions)
    
    [uh,ph,info_stress_next] = vem_mixed_concentration_active_stress_k1_MBC_battery(node,elem,data_stress,bdStruct1,phi0,phi1,info_diffusion);
    
    %% k_2 = 0 for transport mixed formulation (mixed boundary conditions) 
    
    [zetah,phih,info_diffusion_next] = vem_mixed_stress_assisted_diffusion_k1_MBC_battery(node,elem,data_diffusion,bdStruct2,m0,m1,info_stress_next);

    %% Track iterations

    it_ct = it_ct+1;
    it = [it, it_ct];

    %% Stop criteria

    if size(it,2) > 1
        
        %% \norm{sol_{l} - sol_{l-1}} < tol for stop criteria

        u_stop = error_V1_battery(node,elem,uh-uh_prev,info_stress_next,data_stress);
        p_stop = error_Q1_battery(node,elem,ph-ph_prev,data_stress);     
        [zeta_stop,phi_stop] = error_V2_Q2_k1_battery(node,elem,zetah-zetah_prev,phih-phih_prev,info_diffusion_next,data_diffusion,M,m0,m1,info_stress_next);
        sqrt(u_stop.^2+p_stop.^2)+sqrt(zeta_stop.^2+phi_stop.^2)
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

%% Projection for plotting

[zetahI,phihI,nodeI_zeta,elemI_zeta] = projection_V2_k1(node,elem,zetah,phih,info_diffusion_next,data_diffusion,uh);
[phI,nodeI_p,elemI_p] = projection_Q1(node,elem,ph,1,uh);

% Plots

figure(1)
showresult_u_battery(node,elem,uh./1.e-6);

figure(2); 
showresult_zeta_battery(nodeI_zeta,elemI_zeta,zetahI);


figure(3);
showresult_p_phi_battery(nodeI_p,elemI_p,phI,phihI);

drawnow;   

toc;

%% Save data for concentration along radial axis

phih1 = phih; 

% run time check

tic;

%% Neumann condition

bdNeumann1 = '(x.^2+y.^2)==5.^2'; %'x==1 | y==1 | x==0 | y==0';%
bdNeumann2 = '(x.^2+y.^2)==1.^2'; % '(x.^2+y.^2)==1.^2'; %

%% Physical parameters

[lambda,mu,m0,m1,phi0,phi1,theta] = parameters_battery();

%% stress-assisted coefficient

m0 = 1;
m1 = 1.e-5;

%% Data for vem

data_stress = mixed_concentration_active_stress_data_MBC_battery(mu,lambda);
data_diffusion = mixed_stress_assisted_diffusion_data_MBC_battery(theta,mu);

%% Obtain M parameter using optimization routine

M = 1;

fprintf('Value of M: %d\n', M)

%% Stop tolerance for Piccard iteration

tol = 1.e-8;
    
%% Load mesh and domain

mesh_file = 'meshes/circle_circle_10000.mat';
load(mesh_file);

%% print actual mesh and # elements

fprintf(strcat('mesh: ',mesh_file, '\n'));
h = size(elem,1);
fprintf('Number of elements: %d\n', h)

%% Set boundary information

bdStruct1 = setboundary_curve(1/sqrt(sqrt(h)),node,elem,bdNeumann1);
bdStruct2 = setboundary_curve(1/sqrt(sqrt(h)),node,elem,bdNeumann2);

showmesh(node,elem);

%% initial conditions:
%% Coefficients of a polynom in \mathcal{M}_{1}

info_diffusion = [1;0;0];

%% Tracking iterations

it = [];
it_ct = 0;

%% Fixed point iteration

while true

    %% Solve using vem

    %% k_1 = 2 for elasticity mixed formulation (mixed boundary conditions)
    
    [uh,ph,info_stress_next] = vem_mixed_concentration_active_stress_k1_MBC_battery(node,elem,data_stress,bdStruct1,phi0,phi1,info_diffusion);
    
    %% k_2 = 0 for transport mixed formulation (mixed boundary conditions) 
    
    [zetah,phih,info_diffusion_next] = vem_mixed_stress_assisted_diffusion_k1_MBC_battery(node,elem,data_diffusion,bdStruct2,m0,m1,info_stress_next);

    %% Track iterations

    it_ct = it_ct+1;
    it = [it, it_ct];

    %% Stop criteria

    if size(it,2) > 1
        
        %% \norm{sol_{l} - sol_{l-1}} < tol for stop criteria

        u_stop = error_V1_battery(node,elem,uh-uh_prev,info_stress_next,data_stress);
        p_stop = error_Q1_battery(node,elem,ph-ph_prev,data_stress);     
        [zeta_stop,phi_stop] = error_V2_Q2_k1_battery(node,elem,zetah-zetah_prev,phih-phih_prev,info_diffusion_next,data_diffusion,M,m0,m1,info_stress_next);
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

%% Projection for plotting

[zetahI,phihI,nodeI_zeta,elemI_zeta] = projection_V2_k1(node,elem,zetah,phih,info_diffusion_next,data_diffusion,uh);
[phI,nodeI_p,elemI_p] = projection_Q1(node,elem,ph,1,uh);

% Plots

figure(4)
showresult_u_battery(node,elem,uh); 

figure(5); 
showresult_zeta_battery(nodeI_zeta,elemI_zeta,zetahI);

figure(6);
showresult_p_phi_battery(nodeI_p,elemI_p,phI,phihI);

drawnow;   

toc;

%% Save data for concentration along radial axis

phih2 = phih; 

%% Plot concentration in radial direction

% Generate n equally spaced points along the line
x_values = linspace(1, max(nodeI_p(:,1)), 10000);
y_values = linspace(0, 0, 10000);
NT = size(elemI_p,1);
concentration1 = zeros(size(y_values));
concentration2 = zeros(size(y_values));

for point = 1:size(x_values,2)-1
    for iel = 1:NT
        index = elemI_p{iel}; 
        x_vertices = nodeI_p(index,1);
        y_vertices = nodeI_p(index,2);
        is_inside = inpolygon(x_values(point), 0, x_vertices, y_vertices);
        if is_inside
            aux = auxgeometry(nodeI_p,elemI_p);
            node = aux.node; 
            elem = aux.elem;
            centroid = aux.centroid; 
            diameter = aux.diameter; 
            area = aux.area;

            loc_dof = [iel NT+iel 2*NT+iel];
            c_phi1 = [phih1(loc_dof(1),1),phih1(loc_dof(2),1),phih1(loc_dof(3),1)];
            c_phi2 = [phih2(loc_dof(1),1),phih2(loc_dof(2),1),phih2(loc_dof(3),1)];
            xK = centroid(iel,1); 
            yK = centroid(iel,2); 
            hK = diameter(iel);
            
            syms x y;

            m1P1 = 1+0*x;
            m2P1 = (x-xK)./hK; 
            m3P1 = (y-yK)./hK ; 

            phi1 = c_phi1(1)*m1P1+c_phi1(2)*m2P1+c_phi1(3)*m3P1;
        
            phi1 = matlabFunction(phi1,'Vars',{x,y});  

            phi2 = c_phi2(1)*m1P1+c_phi2(2)*m2P1+c_phi2(3)*m3P1;
        
            phi2 = matlabFunction(phi2,'Vars',{x,y});  

            concentration1(point) = phi1(x_values(point),y_values(point));
            concentration2(point) = phi2(x_values(point),y_values(point));
        end
    end
end

figure(242)
hold on
plot(x_values(5:end-5),concentration1(5:end-5),'b-', 'DisplayName', 'Constant stress-assisted coefficient')
plot(x_values(5:end-5),concentration2(5:end-5),'r-', 'DisplayName', 'Stress-assisted coefficient')
xlabel('Radius');
ylabel('Concentration');
title('Plot of Concentration in the radial direction');
legend('show');
hold off