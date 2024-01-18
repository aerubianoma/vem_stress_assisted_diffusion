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

%% For robust analysis

%% For $\lambda$

%variation = [5.e1;1.e2;5.e2;1.e3;5.e3];

%% For $\mu$

%variation = [5.e1;1.e2;5.e2;1.e3;5.e3];

%% For $\theta$

variation = [5.e-4;1.e-4;5.e-3;1.e-3;5.e-2];

errors_V1_Q1 = zeros(maxIt,size(variation,2));
errors_V2_Q2 = zeros(maxIt,size(variation,2));
M_values = zeros(size(variation));

for c = 1:size(variation,1)

%% Variable to check robustness

lambda = variation(c);
variable_string = '$\lambda=$';

fprintf('constant: %d\n', variation(c));

%% Data for vem

data_stress = mixed_concentration_active_stress_data_MBC(mu,lambda,phi0,phi1);
data_diffusion = mixed_stress_assisted_diffusion_data_MBC(mu,lambda,theta,phi0,phi1,m0,m1);

%% Obtain M parameter using optimization routine

xM = fmincon(@(p) -norm(data_diffusion.M_1(p),'fro'),[0,0],[],[],[],[],[0,0],[1,1]);
M = data_diffusion.M_1(xM);
M = abs(M(1,1));
1/M
theta<1/M
end


