function [zetahI,phihI,nodeI,elemI] = projection_V2_k1(node,elem,zetah,phih,info,pde,uh)

% interpolation_stress_assisted_diffusion returns the piecewise basic data structure of 
% L2 projection of uh and L2 projection of ph for mixed VEMs of the stress assisted diffusion problem.
%
% Copyright (C)  Terence Yu. 
% Modified by Rekha Khot and Andres E. Rubiano

auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; 
edge = auxT.edge;
NE = size(edge,1);

u_1 = uh(1:size(node,1),1);
u_2 = uh(size(node,1)+NE+1:size(node,1)+NE+size(node,1),1);
node = node + [u_1,u_2];

%% nodeI, elemI

nodeI = node(horzcat(elem{:}),:);
elemLen = cellfun('length',elem);
elemI = mat2cell(1:sum(elemLen), 1, elemLen)';

%% Get Ph and chi

Ph = info.Ph;

% elementwise global index

index = info.elem2dof; 

% elementwise numerical d.o.f.s

chi = cellfun(@(id) zetah(id), index, 'UniformOutput', false); 

%% Get auxiliary data

aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
diameter = aux.diameter;
centroid = aux.centroid; 

%% uhI %%

NT = size(elem,1); 
zetahI = cell(NT,1);  
phihI = cell(NT,1);

for iel = 1:NT    

    % element information

    index = elem{iel}; 
    Nv = length(index);   
    x = node(index,1); 
    y = node(index,2);
    xK = centroid(iel,1);
    yK = centroid(iel,2);
    hK = diameter(iel);

    %% Basis for $G_1^\nabla(E)$ %%

    m1N1 = @(x,y) [1/hK + 0*x , 0 + 0*x];
    m2N1 = @(x,y) [0 + 0*x , 1/hK + 0*x];
    m3N1 = @(x,y) [2*(x-xK)./hK^2 + 0*x ,0 + 0*x];
    m4N1 = @(x,y) [(y-yK)./hK^2 + 0*x , (x-xK)./hK^2 + 0*x];
    m5N1 = @(x,y) [0 + 0*x , 2*(y-yK)./hK^2 + 0*x];

    %% Basis for $G_1^\oplus(E)$ %%

    m1S1 = @(x,y) [-1*(y-yK)./hK^2 + 0*x , (x-xK)./hK^2 + 0*x];

    %% Basis for $P_1(E)$

    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK ; 

    m = {m1,m2,m3};
    
    % coefficients of projection
    a = Ph{iel}*chi{iel};  

    % L2 projection of uh

    zetahI{iel} = a(1)*m1N1(x,y)+a(2)*m2N1(x,y)+a(3)*m3N1(x,y)+a(4)*m4N1(x,y)+a(5)*m5N1(x,y)+a(6)*m1S1(x,y);

    % L2 projection of phih

    loc_dof = [iel NT+iel 2*NT+iel];

    phihI{iel} = phih(loc_dof(1),1)*m1(x,y)+phih(loc_dof(2),1)*m2(x,y)+phih(loc_dof(3),1)*m3(x,y);

end

zetahI = vertcat(zetahI{:});
phihI = vertcat(phihI{:});