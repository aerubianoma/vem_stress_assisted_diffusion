function showresult_u_p(nodeu,elemu,nodep,elemp,u,uh,p,ph)
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

% Change here to plot 4 graphs instead of 2

tiledlayout(1,2)

%% Plot exact solution u magnitude

% nexttile;
% u_exact = u(node);
% showsolution(nodeu,elemu,vecnorm(u_exact,2,2),'$|\mathbf{u}_{h}|$');


%% Plot numerical solution u magnitude

nexttile;
uh_plot = [uh(1:size(nodeu,1),1),uh(size(nodeu,1)+NE+1:size(nodeu,1)+NE+size(nodeu,1),1)];
showsolution(nodeu,elemu,vecnorm(uh_plot,2,2),'$|\mathbf{u}_{h}|$');

%% Plot exact solution p

%nexttile;
%p_exact = p(node)
%showsolution(nodep,elemp,p_exact(1:size(node,1)),'$\tilde{p}$');

%% Plot numerical solution p

nexttile;
showsolution(nodep,elemp,ph,'$\tilde{p}_h$');



