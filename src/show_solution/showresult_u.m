function showresult_u(node,elem,u,uh,uname,uhname)
%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

clf;  % clear figures

set(gcf,'Units','normal');
%set(gcf,'Position',[0.25,0.25,0.7,0.25]);

auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; 
edge = auxT.edge;
NE = size(edge,1);

%% Plot mesh

%tiledlayout(1,5)

%nexttile;
%showmesh(node,elem); 
%hold on;
%plot(node(:,1),node(:,2),'k.', 'MarkerSize', 4);

tiledlayout(1,2)

%% Plot exact solution first component

%nexttile;
%u_exact = u(node); u_exact = u_exact(:,1);
%showsolution(node,elem,u_exact(1:size(node,1)),'$u_1$');

%% Plot numerical solution

nexttile;
showsolution(node,elem,uh(1:size(node,1),1),'$u_{1,h}$');

%% Plot exact solution second component

%nexttile;
%u_exact = u(node); u_exact = u_exact(:,2);
%showsolution(node,elem,u_exact(1:size(node,1)),'$u_2$');


%% Plot numerical solution

nexttile;
showsolution(node,elem,uh(size(node,1)+NE+1:size(node,1)+NE+size(node,1),1),'$u_{2,h}$');


%nexttile;
%showmesh(node+[uh(1:size(node,1),1),uh(size(node,1)+NE+1:size(node,1)+NE+size(node,1),1)],elem); 
%hold on;
%plot(uh(1:size(node,1),1),uh(size(node,1)+NE+1:size(node,1)+NE+size(node,1),1),'k.', 'MarkerSize', 4);
