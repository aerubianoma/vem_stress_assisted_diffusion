function [zetah,phih,info] = vem_mixed_stress_assisted_diffusion_k1_MBC(node,elem,pde,bdStruct,m0,m1,info_stress)

%% Get auxiliary data for the domain %%

% Auxgeometry

aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; 
diameter = aux.diameter;
area = aux.area;

% Auxstructure

auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; 
edge = auxT.edge;
edge2elem = auxT.edge2elem;

% Number of elementes and edges

%N = size(node,1);
NT = size(elem,1); 
NE = size(edge,1);

%% Polynomial space features %%

% Considering the space $(P_k(E))^2\subseteq V^{h,k}(E)$:
% \dim(V^{h,k}(E)) = (k+1)*Ne + n_{k-1}^{\nabla} + n_{k}^{\oplus}$, Here we are in the case $k = 1$

%k = 1;
%Ndof = (k+1)*NE + 2 + 1;
%NdofEdges = 2*NE;  
%NdofNabla = 2;  
%NdofPlus = 1;

% For the space $Q^{h,k}$:
% $\dim(Q^{h,k}(E))= \frac{(k+1)(k+2)}{2}$

%NdofQ = 3;


%% Dimension for polinomial basis %%

% $\dim((P_1(E))^2)=2*n_1$ with $n_k = \frac{(k+1)(k+2)}{2}$

dimP1_2 = 6;
dimP1 = 3;
dimGN1 = 5;
%dimGN0 = 2;
dimGS1 = 1;

%% Get elementwise signs of basis functions %%

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
    sgnBase{iel} = [sgnelem;sgnelem;ones(3,1)];
end

%% Compute and assemble the linear system %%

NNdofA = 2*NE+(2+1)*NT; 
NNdofB = 3*NT;  
NNdof = NNdofA + NNdofB;

elemLen = cellfun('length',elem); 

nnzA = sum((2*elemLen+2+1).^2); 
nnzB = sum((2*elemLen+2+1).*dimP1);
nnzC = dimP1.^2;

iiA = zeros(nnzA,1); 
jjA = zeros(nnzA,1); 
ssA = zeros(nnzA,1);

iiB = zeros(nnzB,1); 
jjB = zeros(nnzB,1); 
ssB = zeros(nnzB,1);

iiC = zeros(nnzC,1); 
jjC = zeros(nnzC,1); 
ssC = zeros(nnzC,1);

elemb = zeros(NNdofB,1);  

Gb = zeros(NNdofB,1); 

idA = 0; 
idB = 0; 
idC = 0;  
ib = 0;

% matrix for error evaluation

Ph = cell(NT,1); 
elem2dof = cell(NT,1);
div = cell(NT,1);

%% Get element-wise matrices %%

zeta_exact_global = cell(NT,1);

for iel = 1:NT
    
    %% Element information%%
    
    index = elem{iel};
    indexEdge = elem2edge{iel};
    Nv = length(index);
    Ndof = 2*Nv+2+1;
    xK = centroid(iel,1);
    yK = centroid(iel,2);
    hK = diameter(iel);
    x = node(index,1);
    y = node(index,2); 
    sgn = sgnBase{iel}(1:Nv);
    % Loop index for vertices or edges %
    
    v1 = 1:Nv; v2 = [2:Nv,1]; 
    w1 = zeros(1,Nv); w2 = w1;
    for nn = 1:Nv
        if sgn(nn)==-1
            w1(nn) = v2(nn); 
            w2(nn) = v1(nn);
        else
            w1(nn) = v1(nn);
            w2(nn) = v2(nn);
        end
    end
    
    %% First set of dofs %%
    
    % First set of edge points

    xe1 = -1*(x(w1)-x(w2))*(-1/(2*sqrt(3)))+(x(w1)+x(w2))/2;
    ye1 = -1*(y(w1)-y(w2))*(-1/(2*sqrt(3)))+(y(w1)+y(w2))/2; 
    
    % Second set edge of points

    xe2 = -1*(x(w1)-x(w2))*(1/(2*sqrt(3)))+(x(w1)+x(w2))/2;
    ye2 = -1*(y(w1)-y(w2))*(1/(2*sqrt(3)))+(y(w1)+y(w2))/2; 
    
    % Normal vector to the edge

    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; 
    he = sqrt(sum(Ne.^2, 2));

    % Information for the quadrature
    
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    %% Scaled monomials %%
    
    %% basis for $G_1^\nabla(E)$ %%

    m1N1 = @(x,y) [1/hK ; 0];
    m2N1 = @(x,y) [0 ; 1/hK];
    m3N1 = @(x,y) [2*(x-xK)./hK^2 ; 0];
    m4N1 = @(x,y) [(y-yK)./hK^2 ; (x-xK)./hK^2];
    m5N1 = @(x,y) [0 ; 2*(y-yK)./hK^2];

    %% basis for $G_1^\oplus(E)$ %%

    m1S1 = @(x,y) [-1*(y-yK)./hK^2 ; (x-xK)./hK^2];
    
    %% Base for $(P_1(E))^2$ %%

    basisP1_2 = @(x,y) [m1N1(x,y), m2N1(x,y), m3N1(x,y), m4N1(x,y), m5N1(x,y), m1S1(x,y)];

    %% Second set of dofs %%

    m1N0 = @(x,y) m1N1(x,y);
    m2N0 = @(x,y) m2N1(x,y);

    %% Third set of dofs %%

    %m1S1;
    
    %% Basis for $P_2(E)$ %%
    
    m1P2 = @(x,y) 1;                  
    m2P2 = @(x,y) (x-xK)./hK;             
    m3P2 = @(x,y) (y-yK)./hK;             
    m4P2 = @(x,y) (x-xK).^2/hK^2;         
    m5P2 = @(x,y) (x-xK).*(y-yK)./hK^2;   
    m6P2 = @(x,y) (y-yK).^2./hK^2;  

    basisP2 = @(x,y) [m1P2(x,y), m2P2(x,y), m3P2(x,y), m4P2(x,y), m5P2(x,y), m6P2(x,y)];

    %% Basis for $P_1(E)$ %%

    basisP1 = @(x,y) [m1P2(x,y), m2P2(x,y), m3P2(x,y)];
    
    %% Transition matrix %%

    D = zeros(Ndof,dimP1_2); 
    
    % $D_{ij} = dof_i(g_j)$

    % This part covers the dofs on the edges, for this case the dofs are ordered by:
    % - 1:Nv corresponds to [xe1,ye1] points (first points for quadrature on each edge)
    % - Nv+1:2*Nv corresponds to [xe2,ye2] points (second points for quadrature on each edge)
    
    for i = 1:Nv 
        va = basisP1_2(xe1(i),ye1(i)); 
        vb = basisP1_2(xe2(i),ye2(i)); 
        D(i,:) = (1/he(i))*dot(va,repmat(Ne(i,:)',1,dimP1_2));
        D(i+Nv,:) = (1/he(i))*dot(vb,repmat(Ne(i,:)',1,dimP1_2));
    end

    % This part covers the dofs inside E (moments), these are ordered by
    % - 2*Nv+1 and 2*Nv+2 corresponds to the \nabla part
    % - 2*Nv+3 corresponds to the \osum part

    D(2*Nv+1,:) = (1/area(iel))*integralTri(@(x,y) sum(basisP1_2(x,y).*m1N0(x,y)),2,nodeT,elemT);
    D(2*Nv+2,:) = (1/area(iel))*integralTri(@(x,y) sum(basisP1_2(x,y).*m2N0(x,y)),2,nodeT,elemT);
    D(2*Nv+2+1,:) = (1/area(iel))*integralTri(@(x,y) sum(basisP1_2(x,y).*m1S1(x,y)),2,nodeT,elemT); 
    
    %% $L^2(\Omega)$-projection %%

    % $B^\nabla = -\int_{E}\vdiv \phi_i m_{\beta+1} + \int_{\partial E} \phi_i \cdot n_E m_{\beta+1} = B^\nabla_1 + B^\nabla_2$

    %BNabla = zeros(dimGN1,Ndof);
    %BNabla1 = zeros(dimGN1,Ndof);
    BNabla2 = zeros(dimGN1,Ndof);

    % For $B^\nabla_2$ we know that $\phi_i\cdot n_E$ in $\partial E$ is a polynom of degree $1$.
    % We use gaussian quadrature with two points, every edge, $e_{i}$, is the line $\overline{v_{i}v_{i+1}}$, then its integral can be approximated by:
    %\int_{\partial e_{i}} f(x)\, dx = \frac{h_e}{2}(f(q_1)+f(q_2)),
    % where:
    % $q_1 = \frac{-(v_{i+1}-v_{i})}{2\sqrt{3}}+\frac{v_{i+1}+v_{i}}{2}$,
    % $q_2 = \frac{v_{i+1}-v_{i}}{2\sqrt{3}}+\frac{v_{i+1}+v_{i}}{2}$.
    % With the definition of the dofs on the edges we obtain the following:

   for iedge = 1:Nv
        mPoint1_BNabla2 = basisP2(xe1(iedge),ye1(iedge));
        mPoint2_BNabla2 = basisP2(xe2(iedge),ye2(iedge));
        BNabla2(:,iedge) = (he(iedge)/2)*(mPoint1_BNabla2(2:end)');
        BNabla2(:,iedge+Nv) = (he(iedge)/2)*(mPoint2_BNabla2(2:end)');
   end

    % $H = \int_E m_{\sigma}m_{\tao}$

    %H = zeros(dimP1,dimP1);
    H = integralTri(@(x,y) basisP1(x,y)'*basisP1(x,y),2,nodeT,elemT);
    
    % $W1 = -\int_E \phi_i \cdot \nabla m_\tao$
    % Here we use the definition of dofs relationed with $g_0^{\nabla}$

    W1 = zeros(dimP1,Ndof);
    W1(2,2*Nv+1) = -area(iel);
    W1(3,2*Nv+2) = -area(iel);
    
    % $W2 = \int_{\partial E} \phi_i \cdot$
    % Here is the same reasoning as with B_2^{\nabla}
    
    W2 = zeros(dimP1,Ndof);
    for iedge = 1:Nv
        mPoint1_W2 = basisP1(xe1(iedge),ye1(iedge));
        mPoint2_W2 = basisP1(xe2(iedge),ye2(iedge));
        W2(:,iedge) = (he(iedge)/2)*(mPoint1_W2');
        W2(:,iedge+Nv) = (he(iedge)/2)*(mPoint2_W2');
    end
    
    % $H^# = \int_E m_\sigma m_{\beta+1}$

    %HHash = zeros(dimGN1,dimP1);
    
    mP_beta_1 = @(x,y) [m2P2(x,y), m3P2(x,y), m4P2(x,y), m5P2(x,y), m6P2(x,y)];
    HHash = integralTri(@(x,y) mP_beta_1(x,y)'*basisP1(x,y),2,nodeT,elemT);

    BNabla1 = -HHash*H^-1*(W1+W2);

    BNabla = BNabla1 + BNabla2;

    BSum = zeros(dimGS1,Ndof);
    BSum(1,2*Nv+2+1) = area(iel);

    B = [BNabla;BSum];

    G = B*D;
    
    %% Cross checking for correctness of code %%

    % G1 = integralTri(@(x,y) basisP1_2(x,y)'*basisP1_2(x,y),2,nodeT,elemT);
    % norm(G1-G)
    % 
    %% Sign matrix and sign vector %%
    
    sgnelem = sgnBase{iel}; 
    sgnK = sgnelem*sgnelem'; 

    %% Projection %%

    Pis = G\B;
    Pi  = D*Pis;
    I = eye(size(Pi));  

    %% Approximated solution from previous iteration

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
    Ms = integralTri(@(x,y)(M1([x,y])*basisP1_2(x,y))'*basisP1_2(x,y),4,nodeT,elemT);
    Gs = integralTri(@(x,y) M1([x,y]),4,nodeT,elemT);
    AK  = Pis'*Ms*Pis + norm(Gs,'fro')*(I-Pi)'*(I-Pi);
    AK = AK.*sgnK;
    AK = reshape(AK,1,[]); % straighten as row vector for easy assembly 

    %% Stiffness matrix B %%

    BK = (W1+W2);
    BK = BK.*repmat(sgnelem',dimP1,1);
    BK = reshape(BK,1,[]); % straighten as row vector for easy assembly
    
    %% stiffness matrix C %%
    
    CK = pde.theta*H;
    CK = reshape(CK,1,[]);
    
    %% Assembly index for bilinear forms %%

    indexDofV = [elem2edge{iel},NE+elem2edge{iel},2*NE+iel,2*NE+iel+NT,2*NE+iel+2*NT];
    indexDofQ = [iel,iel+NT,iel+2*NT];

    iiA(idA+1:idA+Ndof^2) = reshape(repmat(indexDofV, Ndof,1), [], 1);
    jjA(idA+1:idA+Ndof^2) = repmat(indexDofV(:), Ndof, 1);
    ssA(idA+1:idA+Ndof^2) = AK(:);
    idA = idA + Ndof^2;

    iiB(idB+1:idB+Ndof*dimP1) = reshape(repmat(indexDofV, dimP1,1), [], 1);
    jjB(idB+1:idB+Ndof*dimP1) = repmat(indexDofQ(:), Ndof, 1);
    ssB(idB+1:idB+Ndof*dimP1) = BK(:);
    idB = idB + Ndof*dimP1;

    iiC(idC+1:idC+dimP1^2) = reshape(repmat(indexDofQ, dimP1,1), [], 1);
    jjC(idC+1:idC+dimP1^2) = repmat(indexDofQ(:), dimP1, 1);
    ssC(idC+1:idC+dimP1^2) = CK(:);
    idC = idC + dimP1^2;

    %% load vector g %%

    g = pde.g;
    rhs = integralTri(@(x,y) g([x,y])*basisP1(x,y),4,nodeT,elemT); 
    gK = rhs;

    %% Assembly index for rhs %%

    elemb(ib+1:ib+dimP1) = indexDofQ(:);
    Gb(ib+1:ib+dimP1) = gK(:);
    ib = ib + dimP1;
    
    %% Matrix for error evaluation %%

    sgnPis = repmat(sgnelem',size(Pis,1),1);
    Ph{iel} = sgnPis.*Pis; 
    elem2dof{iel} = indexDofV;
    div{iel} = H^-1*(W1+W2).*repmat(sgnelem',dimP1,1);

end

A = sparse(iiA,jjA,ssA,NNdofA,NNdofA);
B = sparse(iiB,jjB,ssB,NNdofA,NNdofB);
C = sparse(iiC,jjC,ssC,NNdofB,NNdofB);

GB = accumarray(elemb,Gb,[NNdofB 1]);

%% Get block linear system %%
 
kk = sparse(NNdof,NNdof); 
ff = zeros(NNdof,1);

kk(1:NNdofA,1:NNdofA) = A;
kk(1:NNdofA, (1:NNdofB)+NNdofA) = B;
kk((1:NNdofB)+NNdofA, 1:NNdofA) = B';
kk((1:NNdofB)+NNdofA,(1:NNdofB)+NNdofA) = -C;

ff((1:NNdofB)+NNdofA) = -GB;

%% Apply Dirichlet boundary conditions %%

bdEdgeD = bdStruct.bdEdgeD;
bdEdgeIdxD = bdStruct.bdEdgeIdxD;
idD = [bdEdgeIdxD, bdEdgeIdxD+NE];
isBdDofD = false(NNdof,1); 
isBdDofD(idD) = true;
bdDofD = (isBdDofD); 
%freeDofD = (~isBdDofD);


% bdval

phi_D = @(p) pde.phi_exact(p);

bd1D = node(bdEdgeD(:,1),:);
bd2D = node(bdEdgeD(:,2),:);

% First set of edge points

bde1D = -1*(bd1D-bd2D)*(-1/(2*sqrt(3)))+(bd1D+bd2D)/2;

% Second set edge of points

bde2D = -1*(bd1D-bd2D)*(1/(2*sqrt(3)))+(bd1D+bd2D)/2;

% Normal vector to the edge

bdNeD = [bd2D(:,2)-bd1D(:,2), bd1D(:,1)-bd2D(:,1)];
bdheD = sqrt(sum(bdNeD.^2, 2));

% Boundary condition

ff(bdDofD) = [(bdheD/2).*phi_D(bde1D);(bdheD/2).*phi_D(bde2D)];

%% Assemble Neumann boundary conditions

bdEdgeN = bdStruct.bdEdgeN;  
bdEdgeIdxN = bdStruct.bdEdgeIdxN;

if ~isempty(bdEdgeN)
    id = [bdEdgeIdxN, bdEdgeIdxN+NE];
    isBdDofN = false(NNdof,1); 
    isBdDofN(id) = true;
    bdDofN = (isBdDofN); 
    freeDofN = (~isBdDofN);

    zeta_N = @(p) pde.zeta_exact(p);

    bd1N = node(bdEdgeN(:,1),:);
    bd2N = node(bdEdgeN(:,2),:);
    
    % First set of edge points
    
    bde1N = -1*(bd1N-bd2N)*(-1/(2*sqrt(3)))+(bd1N+bd2N)/2;
    
    % Second set edge of points
    
    bde2N = -1*(bd1N-bd2N)*(1/(2*sqrt(3)))+(bd1N+bd2N)/2;
    
    % Normal vector to the edge
    
    bdNeN = [bd2N(:,2)-bd1N(:,2), bd1N(:,1)-bd2N(:,1)];
    bdheN = sqrt(sum(bdNeN.^2, 2));
    
    sol = zeros(NNdof,1);
    sol(bdDofN) = [sum(zeta_N(bde1N).*((1./bdheN).*bdNeN),2);sum(zeta_N(bde2N).*((1./bdheN).*bdNeN),2)];
    
    ff = ff - kk*sol;
end

sol(freeDofN) = kk(freeDofN,freeDofN)\ff(freeDofN);

zetah = sol(1:NNdofA); % u = [u1,u2]
phih = sol(NNdofA+1:end);

%% Store information for computing errors %%

info.Ph = Ph;
info.elem2dof = elem2dof; 
info.kk = kk; 
info.zetah = zetah;
info.phih = phih;
info.div = div;

end
