function [zetah,phih,info] = vem_mixed_stress_assisted_diffusion_k0_MBC(node,elem,pde,bdStruct,m0,m1,info_stress)

%% Get auxiliary data

% auxgeometry

aux = auxgeometry(node,elem);
node = aux.node;
elem = aux.elem;
centroid = aux.centroid; 
diameter = aux.diameter; 
area = aux.area;

% auxstructure

auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; 
edge = auxT.edge;

% number

NT = size(elem,1); 
NE = size(edge,1);
NNdofA = NE;  
NNdofB = NT;
NNdof = NNdofA + NNdofB;
Nm = 3; 
Nmh = Nm-1;

%% Get elementwise signs of basis functions

bdEdgeIdx = bdStruct.bdEdgeIdx;  
E = false(NE,1); 
E(bdEdgeIdx) = 1;
sgnBase = cell(NT,1);

for iel = 1:NT
    index = elem{iel};   
    Nv = length(index);
    NdofA = Nv;
    sgnedge = sign(diff(index([1:Nv,1])));
    id = elem2edge{iel};
    sgnbd = E(id);
    sgnedge(sgnbd) = 1;
    sgnelem = ones(NdofA,1);
    sgnelem(1:Nv) = sgnedge; 
    sgnBase{iel} = sgnelem;
end

%% Compute and assemble the linear system

elemLen = cellfun('length',elem); 
nnzA = sum((elemLen).^2);  nnzB = sum((elemLen)*1); nnzC = 1;
iiA = zeros(nnzA,1); jjA = zeros(nnzA,1); ssA = zeros(nnzA,1);
iiB = zeros(nnzB,1); jjB = zeros(nnzB,1); ssB = zeros(nnzB,1);
iiC = zeros(nnzC,1); jjC = zeros(nnzC,1); ssC = zeros(nnzC,1);
elemb = zeros(NNdofB,1);  Fb = zeros(NNdofB,1); 
idA = 0; idB = 0; idC = 0;  ib = 0;
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
div = cell(NT,1);
d = zeros(NT, 1); % for Lagrange multiplier
drhs = 0;

for iel = 1:NT
    
    % ------- element information --------
    index = elem{iel};    
    Nv = length(index);
    NdofA = Nv;
    xK = centroid(iel,1);
    yK = centroid(iel,2);
    hK = diameter(iel);
    x = node(index,1);
    y = node(index,2); 
    v1 = 1:Nv; 
    v2 = [2:Nv,1]; % loop index for vertices or edges
    %xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    %Te = [-Ne(:,2), Ne(:,1)]; % he*te
    he = norm(Ne);
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials --------     
    
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;
    m = @(x,y) [m2(x,y),m3(x,y)];

    basisP0 = @(x,y) [1 + 0*x];

    % -------- transition matrix ----------
    
    D = zeros(NdofA,Nmh); 
    D = Ne;
    
    % --------- elliptic projection -----------
    
    % First term
    
    B = zeros(Nmh,NdofA);
    pm = integralTri(m,2,nodeT,elemT);
    B = -area(iel)^(-1)*hK*repmat(he',2,1).*repmat(pm',1,Nv);
    
    % Second term
    
    bd = 0.5*(m(x(v1),y(v1))+m(x(v2),y(v2)));
    B = B + hK*bd';
    G = B*D;
    
    %%% One can cross check that G = G1 %%%%%%%%%%
    
    % G1 = area(iel)*eye(2);
    
    % -------- sign matrix and sign vector -------
    
    sgnelem = sgnBase{iel}; 
    sgnK = sgnelem*sgnelem'; 
    
    % ------------- stiffness matrix -------------
    
    % projection
    Pis = G\B;   
    Pi  = D*Pis; 
    I = eye(size(Pi));    

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
        
        chi = cellfun(@(id) uh_previous(id), index_local, 'UniformOutput', false);

        c_w = Ph_local{iel}*chi{iel};

        ph_previous = reshape(ph_previous, [], 3);

        c_r = [ph_previous(iel,1),ph_previous(iel,2),ph_previous(iel,3)];

    end

    %% Stiffness matrix A %%
    
    M1_l = M1_local(m0,m1,pde.mu,c_w,c_r,hK,xK,yK); % Here comes the function to calculate locally M1
    M1 = M1_l.M1;
    Gs = integralTri(@(x,y) M1([x,y]),4,nodeT,elemT);
    AK  = Pis'*Gs*Pis + norm(Gs,'fro')*(I-Pi)'*(I-Pi);
    AK = AK.*sgnK;
    AK = reshape(AK',1,[]); % straighten as row vector for easy assembly 
    
    % stiffness matrix B
    
    BK = zeros(NdofA,1); BK(1:Nv) = 1;
    BK = BK.*sgnelem;
    
    H = integralTri(@(x,y) basisP0(x,y)'*basisP0(x,y),2,nodeT,elemT);
    div{iel} = H^-1*BK;

    BK = reshape(BK',1,[]); % straighten as row vector for easy assembly
    
    % stiffness matrix C
    
    CK = pde.theta*area(iel);
    
    % -------- load vector f ----------
    
    gxy = @(x,y) pde.g([x,y]); % g(p) = g([x,y])
    rhs = integralTri(gxy,3,nodeT,elemT); rhs = rhs';
    gK = rhs;
    
    % ------ assembly index for bilinear forms --------
    
    NdofA = Nv;
    NdofB = 1;
    indexDofA = [elem2edge{iel}];
    indexDofB = iel;
    
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
    
    elemb(ib+1:ib+NdofB) = indexDofB(:);
    Fb(ib+1:ib+NdofB) = gK(:);
    ib = ib + NdofB;

    % ------------ Lagrange multiplier -------------
    
    d(iel,:) = integralTri(@(x,y) 1+0*x,3,nodeT,elemT);
    ph = @(x,y) pde.phi_exact([x,y]);
    drhs = drhs+integralTri(ph,4,nodeT,elemT);
    
    % ------- matrix for error evaluation -------
    
    sgnPis = repmat(sgnelem',size(Pis,1),1);
    Ph{iel} = sgnPis.*Pis; 
    elem2dof{iel} = indexDofA;

end

A = sparse(iiA,jjA,ssA,NNdofA,NNdofA);
B = sparse(iiB,jjB,ssB,NNdofA,NNdofB);
C = sparse(iiC,jjC,ssC,NNdofB,NNdofB);
d = d(:); % for Lagrange multiplier
FB = accumarray(elemb,Fb,[NNdofB 1]);

%% Get block linear system

kk = sparse(NNdof+1,NNdof+1);  ff = zeros(NNdof+1,1);
kk(1:NNdofA,1:NNdofA) = A;
kk(1:NNdofA, (1:NNdofB)+NNdofA) = B;
kk((1:NNdofB)+NNdofA, 1:NNdofA) = B';
kk((1:NNdofB)+NNdofA,(1:NNdofB)+NNdofA) = -C;
kk((1:NNdofB)+NNdofA, end) = d;
kk(end, (1:NNdofB)+NNdofA) = d';
ff((1:NNdofB)+NNdofA) = -FB;
ff(end) = drhs;

%% Apply Dirichlet boundary conditions %%

bdEdgeD = bdStruct.bdEdgeD;
bdEdgeIdxD = bdStruct.bdEdgeIdxD;
idD = [bdEdgeIdxD];
isBdDofD = false(NNdof,1); 
isBdDofD(idD) = true;
bdDofD = (isBdDofD); 
%freeDofD = (~isBdDofD);

% bdval

p = pde.phi_exact;
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:); ze = (z1+z2)./2;
eD = z1-z2;  % e = z2-z1
NeD = [-eD(:,2),eD(:,1)];
ff(bdDofD) = p(ze);

%% Assemble Neumann boundary conditions

bdEdgeN = bdStruct.bdEdgeN;  
bdEdgeIdxN = bdStruct.bdEdgeIdxN;

if ~isempty(bdEdgeN) 
    id = [bdEdgeIdxN];
    isBdDofN = false(NNdof,1); 
    isBdDofN(id) = true;
    bdDofN = (isBdDofN); 
    freeDofN = (~isBdDofN);
    u =  pde.zeta_exact;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:); ze = (z1+z2)./2;
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)];
    bdheN = sqrt(sum(Ne.^2, 2));
    bdval = 1/6*sum((u(z1)+4*u(ze)+u(z2)).*Ne,2); % u*n = g
    % sol
    sol = zeros(NNdof+1,1); 
    sol(bdDofN) = bdval;
    ff = ff - kk*sol;
end

%% Set solver

sol(freeDofN) = kk(freeDofN,freeDofN)\ff(freeDofN);
zetah = sol(1:NNdofA); % u = [u1,u2]
phih = sol(NNdofA+1:end-1);

%% Store information for computing errors

info.Ph = Ph; 
info.elem2dof = elem2dof; 
info.kk = kk;
%info.freeDof = freeDof;
info.div = div;
info.phih = phih;
info.zetah = zetah;