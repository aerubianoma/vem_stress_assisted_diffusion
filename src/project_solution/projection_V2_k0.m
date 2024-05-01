function [uhI,phI,nodeI,elemI] = projection_V2_k0(node,elem,zetah,phih,info,uh)
% ProjectionDarcy returns the piecewise basic data structure of elliptic projection 
% of uh and L2 projection of ph for mixed VEMs of the Darcy problem.
%
% Copyright (C)  Terence Yu. 

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
index = info.elem2dof; % elementwise global index
chi = cellfun(@(id) zetah(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s

if length(phih)>size(elem,1)  % lifting mixed vem
    phih = reshape(phih, [], 3);
end

%% Get auxiliary data

aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
%centroid = aux.centroid;
%diameter = aux.diameter; 

%% uhI

NT = size(elem,1); %K = pde.K;  
uhI = cell(NT,1);
phI = cell(NT,1);

for iel = 1:NT    

    % element information
    
    index = elem{iel};  Nv = length(index);   
    %xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    
    % \hat{m}_a = K*grad(hK*m_{a+1})
    
    mh1 = @(x,y) [1+0*x,0+0*x];
    mh2 = @(x,y) [0+0*x,1+0*x];
    
    % coefficients of elliptic projection
    
    a = Ph{iel}*chi{iel};    
    
    % elliptic projection of uh
    
    uhI{iel} = a(1)*mh1(x,y)+a(2)*mh2(x,y);
    
    % L2 projection of ph
    
    phI{iel} = phih(iel,1)*ones(Nv,1);

end

uhI = vertcat(uhI{:});  % uh = [u1,u2]
phI = vertcat(phI{:});