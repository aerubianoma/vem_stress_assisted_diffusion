function showresult_u_p_battery(nodeu,elemu,nodep,elemp,uh,ph)
%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

clf;  % clear figures

set(gcf,'Units','normal');
%set(gcf,'Position',[0.25,0.25,0.7,0.25]);

auxT = auxstructure(nodeu,elemu);
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


%% Plot numerical solution

nexttile;
uh_plot = [uh(1:size(nodeu,1),1),uh(size(nodeu,1)+NE+1:size(nodeu,1)+NE+size(nodeu,1),1)];
showsolution(nodeu,elemu,vecnorm(uh_plot,2,2),'$|\mathbf{u}_{h}|$');


%% Plot numerical solution p

nexttile;
showsolution(nodep,elemp,ph,'$\tilde{p}_h$');
