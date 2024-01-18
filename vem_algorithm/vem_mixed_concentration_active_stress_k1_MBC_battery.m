function [uh,ph,info] = vem_mixed_concentration_active_stress_k1_MBC_battery(node,elem,pde,bdStruct,phi0,phi1,info_diffusion)

%Stokes_mixedVEM_New solves linear Stokes-elasticity using the mixed virtual element method
% in the lowest order k=2:  u = [u1,u2], f = [f1,f2]
%     
%     - mu*div(\eps(u)) + \nabla(p) = f   in Omega
%     p = -div(u) + g  in Omega
%     Dirichlet boundary condition   u = uD on \Gamma
%
%  The mixed formulation is: Find (u,p) \in (V,Q) such that
%
%      a(u,v) + b(v,p) = (f,v),   v in V,
%      b(u,q) + c(p,q) = -(g,q), q in Q,
%
% where,
% - V = (H_0^1(Omega))^2,   Q = L^2(Omega)
% - u is approximated by the lowest order virtual element (k = 2)
% - p by piecewise linear element (k = 2)
% - Q_{k-1}(K) = P_{k-1}(K).
% Copyright (C)  Rekha Khot (modified Stokes VEM code of Terence Yu)

%% Get auxiliary data

% Auxgeometry

aux = auxgeometry(node,elem);
node = aux.node; 
elem = aux.elem;
centroid = aux.centroid; 
diameter = aux.diameter; 
area = aux.area;

% Auxstructure

auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge;
edge = auxT.edge;

% number

N = size(node,1); 
NT = size(elem,1); 
NE = size(edge,1);
NNdofA = (N+NE)*2 + 2*NT;  
NNdofB = 3*NT;  
NNdof = NNdofA + NNdofB;
Nm = 6; Nmm = 2*Nm;
ff = zeros(NNdof,1);

%% Compute and assemble the linear system

elemLen = cellfun('length',elem); 
nnzA = sum((4*elemLen+2).^2);  nnzB = sum((4*elemLen+2)*3); nnzC = 9*NT;
iiA = zeros(nnzA,1); jjA = zeros(nnzA,1); ssA = zeros(nnzA,1);
iiB = zeros(nnzB,1); jjB = zeros(nnzB,1); ssB = zeros(nnzB,1);
iiC = zeros(nnzC,1); jjC = zeros(nnzC,1); ssC = zeros(nnzC,1);
elemb = zeros(NNdofA,1);  Fb = zeros(NNdofA,1); 
elemc = zeros(NNdofB,1);  Fc = zeros(NNdofB,1);
idA = 0; idB = 0; idC = 0; ib = 0; ic = 0;

% matrix for error evaluation

Ph = cell(NT,1); 
elem2dof = cell(NT,1);


for iel = 1:NT

    % ------- element information ----------

    index = elem{iel};     
    Nv = length(index);
    xK = centroid(iel,1); 
    yK = centroid(iel,2); 
    hK = diameter(iel);
    x = node(index,1); 
    y = node(index,2);

    % loop index for vertices or edges

    v1 = 1:Nv; v2 = [2:Nv,1]; 

    % mid-edge points

    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; 

    % he*ne

    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; 

    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % --------------- scaled monomials -----------------

    m1 = @(x,y) 1+0*x;                  gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;             gradm2 = @(x,y) [1/hK+0*x, 0+0*x];
    m3 = @(x,y) (y-yK)./hK;             gradm3 = @(x,y) [0+0*x, 1/hK+0*x];
    m4 = @(x,y) (x-xK).^2/hK^2;         gradm4 = @(x,y) [2*(x-xK)./hK^2, 0+0*x];
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2;   gradm5 = @(x,y) [(y-yK)./hK^2, (x-xK)./hK^2];
    m6 = @(x,y) (y-yK).^2./hK^2;        gradm6 = @(x,y) [0+0*x, 2*(y-yK)./hK^2];

    e11 = @(x,y) [0+0*x, 0+0*x];                    e21 = @(x,y) [0+0*x, 0+0*x];
    e12 = @(x,y) [1/hK+0*x, 0+0*x];                 e22 = @(x,y) [0+0*x, 0+0*x];
    e13 = @(x,y) 0.5*[0+0*x, 1/hK+0*x];             e23 = @(x,y) 0.5*[1/hK+0*x,0+0*x];
    e14 = @(x,y) [2*(x-xK)./hK^2, 0+0*x];           e24 = @(x,y) [0+0*x, 0+0*x];
    e15 = @(x,y) [(y-yK)./hK^2, (x-xK)./(2*hK^2)];  e25 = @(x,y) [(x-xK)./(2*hK^2),0+0*x];
    e16 = @(x,y) [0+0*x, (y-yK)./hK^2];             e26 = @(x,y) [(y-yK)./hK^2,0+0*x];
    e17 = @(x,y) [0+0*x, 0+0*x];                    e27 = @(x,y) [0+0*x, 0+0*x];
    e18 = @(x,y) 0.5*[0+0*x, 1/hK+0*x];             e28 = @(x,y) 0.5*[1/hK+0*x,0+0*x];
    e19 = @(x,y) [0+0*x, 0+0*x];                    e29 = @(x,y) [0+0*x, 1/hK+0*x];
    e110 = @(x,y) [0+0*x,(x-xK)./hK^2];             e210 = @(x,y) [(x-xK)./hK^2,0+0*x];
    e111 = @(x,y) [0+0*x,(y-yK)./(2*hK^2)];         e211 = @(x,y) [(y-yK)./(2*hK^2),(x-xK)./hK^2];
    e112 = @(x,y) [0+0*x, 0+0*x];                   e212 = @(x,y) [0+0*x, 2*(y-yK)./hK^2];
    
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)];
    Gradm = {gradm1, gradm2, gradm3, gradm4, gradm5, gradm6};    
    divmm = @(x,y) [0+0*x, 1/hK+0*x, 0+0*x, 2*(x-xK)./hK^2, (y-yK)./hK^2, 0+0*x, 0+0*x, 0+0*x, 1/hK+0*x, 0+0*x, (x-xK)./hK^2, 2*(y-yK)./hK^2];
    
    % -------- transition matrix D----------

    NdofBd = 2*Nv; NdofA = 2*NdofBd+2;
    divmm2 = @(x,y) divmm(x,y).*repmat(m2(x,y),1,Nmm);
    divmm3 = @(x,y) divmm(x,y).*repmat(m3(x,y),1,Nmm);
    D = zeros(NdofA, Nmm);
    Dbd = [m(x,y); m(xe,ye)];
    D(1:4*Nv, :) = blkdiag(Dbd, Dbd);
    D(end-1,:) = integralTri(divmm2,4,nodeT,elemT);
    D(end,:) = integralTri(divmm3,4,nodeT,elemT);
    
    % --------- elliptic projection -----------

    % B

    % --- first term ---

    Lapm1 = [0, 0, 0, 2/hK^2, 0, 1/hK^2];
    Lapm2 = [0, 0, 0, 0, 1/(2*hK^2), 0];
    Lapm3 = [0, 0, 0, 1/hK^2, 0, 2/hK^2];
    B = zeros(Nmm, NdofA);
    B(1:Nm, end-1) = 2*pde.mu*hK*Lapm1;
    B(1:Nm, end) = 2*pde.mu*hK*Lapm2;
    B(Nm+1:end, end) = 2*pde.mu*hK*Lapm3;
    B(Nm+1:end, end-1) = 2*pde.mu*hK*Lapm2;

    % --- second term ---

    elem1 = [v1(:), v2(:), v1(:)+Nv]; % elem2dof for [ae, be, me]
    Gradmm = cell(2,Nmm); Emm = cell(2,Nmm);
    Gradmm(1,1:Nm) = Gradm; 
    Gradmm(2,Nm+1:end) = Gradm;
    Emm(1,:) = {e11,e12,e13,e14,e15,e16,e17,e18,e19,e110,e111,e112};
    Emm(2,:) = {e21,e22,e23,e24,e25,e26,e27,e28,e29,e210,e211,e212};

    for im = 1:Nm
        Gradmm{1,im+Nm} = @(x,y) [0+0*x, 0+0*x];
        Gradmm{2,im} = @(x,y) [0+0*x, 0+0*x];
    end

    qmm = cell(1,Nmm);  % q = mu*hK(c2*m2 + c3*m3)
    c2 = [Lapm1, Lapm2];
    c3 = [Lapm2, Lapm3];

    for im = 1:Nmm
        qmm{im} = @(x,y) 2*pde.mu*hK*(c2(im)*m2(x,y) + c3(im)*m3(x,y));
    end

    for s = 1:2
        id = (1:NdofBd) + (s-1)*NdofBd;
        for im = 1:Nmm
            pm = @(x,y) 2*pde.mu*Emm{s,im}(x,y);
            qa = @(x,y) qmm{im}(x,y);
            F1 = 1/6*(sum(pm(x(v1),y(v1)).*Ne, 2) - qa(x(v1),y(v1)).*Ne(:,s));
            F2 = 1/6*(sum(pm(x(v2),y(v2)).*Ne, 2) - qa(x(v2),y(v2)).*Ne(:,s));
            F3 = 4/6*(sum(pm(xe,ye).*Ne, 2) - qa(xe,ye).*Ne(:,s));            
            B(im, id) = accumarray(elem1(:), [F1; F2; F3], [NdofBd, 1]);
        end
    end

    % constraint Bs

    Bs = B;  Bs([1,7,8], :) = 0;
    Bs(1,1:Nv) = ones(1,Nv);
    Bs(7,2*Nv+(1:Nv)) = ones(1,Nv);
    Bs(8,1:Nv) = -y';
    Bs(8,2*Nv+(1:Nv)) = x';

    %  G, Gs

    G = B*D;  Gs = Bs*D;

    %% One can cross-check G = G1 %%

    % G1 = zeros(12,12);
    % for i=1:12
    %    for j = 1:12
    %        pr1 = @(x,y) Emm{1,i}(x,y).*Emm{1,j}(x,y);
    %        pr2 = @(x,y) Emm{2,i}(x,y).*Emm{2,j}(x,y);
    %        I1 = integralTri(pr1,2,nodeT,elemT);
    %        I2 = integralTri(pr2,2,nodeT,elemT);
    %        G1(i,j) = G1(i,j)+sum(I1)+sum(I2);
    %    end
    % end
    % M = (1/area(iel))*integralTri(m,2,nodeT,elemT);
    % G1(1,1:6) = M;
    % G1(7,7:12) = M;

    % ------------- stiffness matrix -------------

    % Projection

    Pis = Gs\Bs;   
    Pi  = D*Pis;   
    I = eye(size(Pi));   

    % Stiffness matrix A for Vh

    AK  = Pis'*G*Pis + (2*pde.mu)*(I-Pi)'*(I-Pi);
    AK = reshape(AK',1,[]); % straighten as row vector for easy assembly

    % Stiffness matrix B for Vh and Qh

    BK = zeros(NdofA,3);
    BK(end-1, 2) = 1;  BK(end, 3) = 1;
    F = 1/6*[(1*Ne); (1*Ne); (4*Ne)]; % [n1, n2]
    BK(1:end-2,1) = accumarray([elem1(:);elem1(:)+NdofBd], F(:), [2*NdofBd 1]);
    BK = reshape(BK',1,[]); % straighten as row vector for easy assembly 

    % Stiffness matrix C for Qh

    CK = zeros(3,3);
    cxy = @(x,y) [m1(x,y);m2(x,y);m3(x,y)]*[m1(x,y),m2(x,y),m3(x,y)];
    CK = (1/pde.lambda)*integralTri(cxy,2,nodeT,elemT);

    % Load vector f and g

    P0K = zeros(2,NdofA);
    P0K(1, end-1) = -1; P0K(2, end) = -1;
    m23 = {m2, m3};

    for s = 1:2
        mc = m23{s};
        F1 = 1/6*(mc(x(v1),y(v1)).*Ne);  % [n1, n2]
        F2 = 1/6*(mc(x(v2),y(v2)).*Ne);
        F3 = 4/6*(mc(xe,ye).*Ne);
        F = [F1; F2; F3];
        P0K(s, 1:NdofBd) = accumarray(elem1(:), F(:,1), [NdofBd 1]);
        P0K(s, NdofBd+1:2*NdofBd) = accumarray(elem1(:), F(:,2), [NdofBd 1]);
    end

    P0K = 1/area(iel)*hK*P0K;
    fxy = @(x,y) pde.f([x,y]); % f(p) = f([x,y]), f = [f1, f2]
    Pf = integralTri(fxy,3,nodeT,elemT); 
    fK = Pf(1)*P0K(1,:) + Pf(2)*P0K(2,:);

    if  ~isa(info_diffusion,'struct')

        c_phi = info_diffusion;

    else
        
        phi_previous = info_diffusion.phih;
        
        % elementwise numerical d.o.f.s
        
        loc_dof = [iel NT+iel 2*NT+iel];

        c_phi = [phi_previous(loc_dof(1),1),phi_previous(loc_dof(2),1),phi_previous(loc_dof(3),1)];

    end
    
    
    l_l = l_local_k_1(phi0,phi1,c_phi,hK,xK,yK);
    l = l_l.g;
    gxy = @(x,y) l([x,y])*[m1(x,y);m2(x,y);m3(x,y)];
    gK = (1/pde.lambda)*integralTri(gxy,4,nodeT,elemT);

    % ------ assembly index for bilinear forms --------

    NdofA = 4*Nv+2; 
    NdofB = 3;

    indexDofA = [elem{iel}, elem2edge{iel}+N, elem{iel}+N+NE, elem2edge{iel}+2*N+NE, iel+2*N+2*NE, iel+2*N+2*NE+NT];
    indexDofB = [iel, iel+NT, iel+2*NT];
    iiA(idA+1:idA+NdofA^2) = reshape(repmat(indexDofA, NdofA,1), [], 1);
    jjA(idA+1:idA+NdofA^2) = repmat(indexDofA(:), NdofA, 1);
    ssA(idA+1:idA+NdofA^2) = AK(:);
    idA = idA + NdofA^2;

    iiB(idB+1:idB+NdofA*NdofB) = reshape(repmat(indexDofA, NdofB,1), [], 1);
    jjB(idB+1:idB+NdofA*NdofB) = repmat(indexDofB(:), NdofA, 1);
    ssB(idB+1:idB+NdofA*NdofB) = BK(:);
    idB = idB + NdofA*NdofB;

    iiC(idC+1:idC+NdofB^2) = reshape(repmat(indexDofB, NdofB,1), [], 1);
    jjC(idC+1:idC+NdofB^2) = repmat(indexDofB(:), NdofB, 1);
    ssC(idC+1:idC+NdofB^2) = CK(:);
    idC = idC + NdofB^2;

    % ------- assembly index for rhs -------

    elemb(ib+1:ib+NdofA) = indexDofA(:);
    Fb(ib+1:ib+NdofA) = fK(:);
    ib = ib + NdofA;

    elemc(ic+1:ic+NdofB) = indexDofB(:);
    Fc(ic+1:ic+NdofB) = gK(:);
    ic = ic + NdofB;
    
    % ------- matrix for error evaluation -------

    Ph{iel} = Pis; 
    elem2dof{iel} = indexDofA;
end

A = sparse(iiA,jjA,ssA,NNdofA,NNdofA);
B = sparse(iiB,jjB,ssB,NNdofA,NNdofB);
C = sparse(iiC,jjC,ssC,NNdofB,NNdofB);
Fb = accumarray(elemb,Fb,[NNdofA 1]);
Fc = accumarray(elemc,Fc,[NNdofB 1]);

%% Get block linear system

kk = sparse(NNdof,NNdof); 

kk(1:NNdofA,1:NNdofA) = A;
kk(1:NNdofA, (1:NNdofB)+NNdofA) = -B;
kk((1:NNdofB)+NNdofA, 1:NNdofA) = -B';
kk((1:NNdofB)+NNdofA, (1:NNdofB)+NNdofA) = -C;

ff(1:NNdofA) = Fb;
ff(NNdofA+(1:NNdofB)) = -Fc;

%% Assemble Neumann boundary conditions

bdEdgeN = bdStruct.bdEdgeN;
bdEdgeIdxN = bdStruct.bdEdgeIdxN;
bdDofN = [bdEdgeN, N+bdEdgeIdxN, N+NE+bdEdgeN, 2*N+NE+bdEdgeIdxN];

if ~isempty(bdEdgeN)
    %g_N = pde.Eu; pe = pde.pexact;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:); ze = 0.5*(z1+z2);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    %pm1 = pde.mu*g_N(z1); pm2 = pde.mu*g_N(z2); pm3 = pde.mu*g_N(ze);
    sigma_N = pde.sigma_N;
    FN1 = -sigma_N(z1).*Ne(:,1); 
    FN2 = -sigma_N(z2).*Ne(:,1);
    FN3 = -sigma_N(ze).*Ne(:,1);
    FN4 = -sigma_N(z1).*Ne(:,1);  
    FN5 = -sigma_N(z2).*Ne(:,1);
    FN6 = -sigma_N(ze).*Ne(:,1);
    FN = [FN1; FN2; FN3; FN4; FN5; FN6];
    ff = ff + accumarray(bdDofN(:), FN(:), [NNdof 1]);
end

%% Apply Dirichlet boundary conditions %%

% bdDof, freeDof

bdNodeIdxD = bdStruct.bdNodeIdxD;bdEdgeIdxD = bdStruct.bdEdgeIdxD;
bdEdgeD = bdStruct.bdEdgeD; 
idD = [bdNodeIdxD; bdEdgeIdxD+N; bdNodeIdxD+N+NE; bdEdgeIdxD+2*N+NE];
isBdDof = false(NNdof,1); isBdDof(idD) = true;
bdDofD = (isBdDof); 
freeDof = (~isBdDof);

% bdval

u_D = pde.u_D;
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:); ze = (z1+z2)./2;
uv = u_D(node(bdNodeIdxD,:));
ue = u_D(ze);
bdval = [uv(:,1); ue(:,1); uv(:,2); ue(:,2)]; 

% sol

sol = zeros(NNdof,1); 
sol(bdDofD) = bdval;

ff = ff - kk*sol;

%% Set solver %%

sol(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
uh = sol(1:NNdofA); % u = [u1(zv); u1(ze); u2(zv); u2(ze); divu*m2; divu*m3]
ph = sol(NNdofA+1:end);

%% Store information for computing errors

info.Ph = Ph; info.elem2dof = elem2dof; 
% info.kk = kk; %info.freeDof = freeDof;
info.D = D;
info.uh = uh;
info.ph = ph;