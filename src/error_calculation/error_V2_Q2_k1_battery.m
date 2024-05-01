function [uL2,pL2] = error_V2_Q2_k1_battery(node,elem,zetah,phih,info,pde,M,m0,m1,info_stress)

% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s 

% coefficient matrix

%M = pde.M; 
%M1 = pde.M_1;
%select = @(M,r,c) M(r,c);

%% Get Ph and chi %%

Ph = info.Ph;

% elemenwise global index

index = info.elem2dof; 

% elementwise numerical d.o.f.s

chi = cellfun(@(id) zetah(id), index, 'UniformOutput', false);

%% Get auxiliary data %%

% exact solution

%zetae = pde.zeta_exact;  phie = pde.phi_exact;
div = info.div;

% auxiliary data structure

aux = auxgeometry(node,elem); elem = aux.elem;
NT = size(elem,1);

%% Compute L2 error %%

ErruL2 = 0; ErrpL2 = 0;
uL2 = 0; pL2 = 0;

for iel = 1:NT

    % element information

    index = elem{iel};  Nv = length(index);    
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); hK = aux.diameter(iel); 

    %% Basis for $G_1^\nabla(E)$ %%

    m1N1 = @(x,y) [1/hK + 0*x , 0 + 0*x];
    m2N1 = @(x,y) [0 + 0*x , 1/hK + 0*x];
    m3N1 = @(x,y) [2*(x-xK)./hK^2 + 0*x ,0 + 0*x];
    m4N1 = @(x,y) [(y-yK)./hK^2 + 0*x , (x-xK)./hK^2 + 0*x];
    m5N1 = @(x,y) [0 + 0*x , 2*(y-yK)./hK^2 + 0*x];

    %% Basis for $G_1^\oplus(E)$ %%

    m1S1 = @(x,y) [-1*(y-yK)./hK^2 + 0*x , (x-xK)./hK^2 + 0*x];
    
    %% Basis for $P_1(E)$

    m1_1 = @(x,y) 1+0*x;
    m2_1 = @(x,y) (x-xK)./hK; 
    m3_1 = @(x,y) (y-yK)./hK;

    %% Divergence basis
    
    %% Basis for $G_1^\nabla(E)$ %%

    dm1N1 = @(x,y) 0;
    dm2N1 = @(x,y) 0;
    dm3N1 = @(x,y) 2./hK^2;
    dm4N1 = @(x,y) 0;
    dm5N1 = @(x,y) 2./hK^2;

    %% Basis for $G_1^\oplus(E)$ %%

    dm1S1 = @(x,y) 0;

    % Projection matrix

    Pis = Ph{iel};

    % Coefficient

    a = Pis*chi{iel};

    % Indexing for phi

    loc_dof = [iel NT+iel 2*NT+iel];

    % Approximation of zeta and phi 

    amK = @(x,y) a(1)*m1N1(x,y)+a(2)*m2N1(x,y)+a(3)*m3N1(x,y)+a(4)*m4N1(x,y)+a(5)*m5N1(x,y)+a(6)*m1S1(x,y);
    pmK = @(x,y) phih(loc_dof(1),1)*m1_1(x,y)+phih(loc_dof(2),1)*m2_1(x,y)+phih(loc_dof(3),1)*m3_1(x,y);
    damk = @(x,y) [m1_1(x,y),m2_1(x,y),m3_1(x,y)]*div{iel}*chi{iel};

    % Errors

    if  ~isa(info_stress,'struct')

        c_w = zeros(1,12);
        c_r = zeros(1,3);

    else

        Ph_local = info_stress.Ph;
        uh_previous = info_stress.uh;
        ph_previous = info_stress.ph;
        
        % elementwise global index
        
        index_local = info_stress.elem2dof; 
        
        % elementwise numerical d.o.f.s
        
        chi2 = cellfun(@(id) uh_previous(id), index_local, 'UniformOutput', false);

        c_w = Ph_local{iel}*chi2{iel};

        ph_previous = reshape(ph_previous, [], 3);

        c_r = [ph_previous(iel,1),ph_previous(iel,2),ph_previous(iel,3)];

    end


    %% Stiffness matrix A %%
    M1_l = M1_local_battery(m0,m1,pde.mu,c_w,c_r,hK,xK,yK); % Here comes the function to calculate locally M1
    M1 = M1_l.M1;

    %errdiv = @(x,y) (pde.div_z([x,y])-damk(x,y)).^2;  
    %errf = @(x,y) (M1([x,y])*(zetae([x,y])-amK(x,y))')'*(zetae([x,y])-amK(x,y))';
    %errp = @(x,y) (phie([x,y])-pmK(x,y)).^2;

    %errdiv = @(x,y) div(x,y).^2;
    
    u = @(x,y) (M1([x,y])*amK(x,y)')'*amK(x,y)';  
    divu = @(x,y) damk(x,y).^2;
    p = @(x,y) pmK(x,y).^2;
    
    % elementwise error

    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    %ErruL2 = ErruL2 + sum(integralTri(errf,4,nodeT,elemT)) + M*integralTri(errdiv,4,nodeT,elemT);
    %ErrpL2 = ErrpL2 + (1/M)*integralTri(errp,4,nodeT,elemT) + pde.theta*integralTri(errp,4,nodeT,elemT);
    uL2 = uL2 + sum(integralTri(u,4,nodeT,elemT)) + M*integralTri(divu,4,nodeT,elemT);
    pL2 = pL2 + (1/M)*integralTri(p,4,nodeT,elemT) + pde.theta*integralTri(p,4,nodeT,elemT);
    
end

h = mean(aux.diameter);
%ErruL2 = sqrt(ErruL2);%/sqrt(uL2);  
%ErrpL2 = sqrt(ErrpL2);%/sqrt(pL2);
uL2 = sqrt(uL2);
pL2 = sqrt(pL2);