function showresult_zeta_battery(node,elem,zetah)
%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

clf;  % clear figures

set(gcf,'Units','normal');
%set(gcf,'Position',[0.25,0.25,0.7,0.25]);

%% Plot mesh

%tiledlayout(1,5)

%nexttile;
%showmesh(node,elem); 
%hold on;
%plot(node(:,1),node(:,2),'k.', 'MarkerSize', 4);

tiledlayout(1,2)

%% Plot exact solution first component

%nexttile;
%zeta_exact = zeta(node); zeta_exact = zeta_exact(:,1);
%showsolution(node,elem,zeta_exact(1:size(node,1)),'$\zeta_1$');


%% Plot numerical solution

nexttile;
showsolution_battery(node,elem,zetah(1:size(node,1),1),'$\zeta_{1,h}$');
hold on 
fplot(@(t) 5*sin(t), @(t) 5*cos(t),'--k','LineWidth',3);
fplot(@(t) sin(t), @(t) cos(t),'--k','LineWidth',3);


%% Plot exact solution second component

%nexttile;
%zeta_exact = zeta(node); zeta_exact = zeta_exact(:,2);
%showsolution(node,elem,zeta_exact(1:size(node,1)),'$\zeta_2$');


%% Plot numerical solution

nexttile;
showsolution_battery(node,elem,zetah(1:size(node,1),2),'$\zeta_{2,h}$');
hold on 
fplot(@(t) 5*sin(t), @(t) 5*cos(t),'--k','LineWidth',3);
fplot(@(t) sin(t), @(t) cos(t),'--k','LineWidth',3);